from math import fmod
from statistics import mean
from constants import SECONDS_IN_A_DAY, DAYS_IN_A_YEAR
import numpy as np
import argparse
from parfile import ParFile
from astropy.time import Time
import sys
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, solar_system_ephemeris, solar_system, AltAz
import astropy.constants as const
from astroplan import Observer, FixedTarget


def bary2topo(source, bary_time):
    lt = bary_time.tdb.light_travel_time(source, location=EarthLocation.of_site("MeerKAT"))
    return (bary_time - lt).utc



def convAngle(str):
    parts=str.split(':')
    if len(parts)!=3:
        raise RuntimeError("unsupported coordinate string {0}".format(str))
    if parts[0].find('-')!=-1:
        return float(parts[0])-float(parts[1])/60.-float(parts[2])/3600.
    else:
        return float(parts[0])+float(parts[1])/60.+float(parts[2])/3600.


def true_anomaly(E, ecc):
    # Return the True Anomaly (in radians) given the Eccentric anomaly
    # (E in radians) and the eccentricity (ecc)
    return 2.0*np.arctan(np.sqrt((1.0+ecc)/(1.0-ecc))*np.tan(E/2.0))

def bisect(func, par_file, lox, hix, TOL=1e-14, MAXIT=200):
    """
    bisect(func, lox, hix, TOL=1e-14, MAXIT=200):
       Try to find a root between 'lox' and 'hix' using a simple
       bisection of the region.  'TOL' is an _absolute_
       tolerance.  'MAXIT' is the maximum number of iterations
    """
    f = func(lox, par_file)
    fmid = func(hix, par_file)
    if (f * fmid >= 0.0):
        print("Root must be bracketed in bisect()!")
        return 0.0
    if (f < 0.0):
        dx, rtb = hix - lox, lox
    else:
        dx, rtb = lox - hix, hix
    for i in range(MAXIT):
        dx = dx * 0.5
        xmid = rtb + dx
        fmid = func(xmid, par_file)
        if (fmid <= 0.0):
            rtb = xmid
        if (abs(dx) < TOL or fmid == 0.0):
            return rtb
    print("Too many bisections in bisect()!")
    return 0.0


def solve_T_at_ntrue_anom(T, par_file, nu=np.pi/2.0):
    return calc_true_anom(T, par_file) - nu

def get_markley_starter(mean_anomaly, ecc, ome):
    FACTOR1 = 3 * np.pi / ( np.pi - 6 /  np.pi)
    FACTOR2 = 1.6 / ( np.pi - 6 / np.pi)

    msq = mean_anomaly * mean_anomaly
    alpha =  FACTOR1 +  FACTOR2 * (np.pi - mean_anomaly) / (1 + ecc)
    d = 3 * ome + alpha * ecc
    alphad = alpha * d
    r = (3 * alphad * (d - ome) + msq) * mean_anomaly
    q = 2 * alphad * ome - msq
    qsq = q * q
    w = pow(np.abs(r) + np.sqrt(qsq * q + r * r), 2.0 / 3)
    return (2 * r * w / (w * w + w * q + qsq) + mean_anomaly) / d


def refine_estimate(mean_anomaly, ecc, ome, E):
    sE = E - np.sin(E)
    cE = 1 - np.cos(E)
    f_0 = ecc * sE + E * ome - mean_anomaly
    f_1 = ecc * cE + ome
    f_2 = ecc * (E - sE)
    f_3 = 1 - f_1
    d_3 = -f_0 / (f_1 - 0.5 * f_0 * f_2 / f_1)
    d_4 = -f_0 / (f_1 + 0.5 * d_3 * f_2 + (d_3 * d_3) * f_3 / 6)
    d_42 = d_4 * d_4
    dE = -f_0 / (f_1 + 0.5 * d_4 * f_2 + d_4 * d_4 * f_3 / 6 - d_42 * d_4 * f_2 / 24)
    return E + dE


def keplers_eqn(mean_anomaly, ecc, tolerance):
    M = fmod(mean_anomaly, 2 * np.pi)
    ome = 1 - ecc
    E = get_markley_starter(mean_anomaly, ecc, ome)
    iter=0
    while True:
        newE = refine_estimate(mean_anomaly, ecc, ome, E)
        iter = iter + 1
        if (E - newE) < tolerance or iter > 1000:
            break
        E = newE

    return E



def calc_true_anom(epoch, par_file):
    difft = (epoch-par_file.get_parameter('T0').value)*SECONDS_IN_A_DAY
    om = par_file.get_parameter('OM').value 
    if(par_file.has_parameter('OMDOT')):
         om = om +  par_file.get_parameter('OMDOT').value * difft/(SECONDS_IN_A_DAY * DAYS_IN_A_YEAR)
    om *= np.pi/180.0
    om = np.fmod(om, np.pi*2.0)
    if (om < 0.0):
        om += (np.pi*2.0)
    time_since_peri = np.fmod(difft, par_file.get_parameter('PB').value * SECONDS_IN_A_DAY)
    if (time_since_peri < 0.0):
        time_since_peri += (par_file.get_parameter('PB').value * SECONDS_IN_A_DAY)
    E = keplers_eqn(time_since_peri * par_file.get_parameter('N_ORBIT').value  , par_file.get_parameter('ECC').value, 1.0E-9)
    ret = np.fmod(true_anomaly(E, par_file.get_parameter('ECC').value) + om, 2*np.pi)
    if (ret < 0.0):
        ret += (np.pi*2.0)
    #print("True anomaly (deg) at ", epoch, " = ", (ret)*180.0/np.pi)

    return ret


def parse_args():
    parser = argparse.ArgumentParser(
        description='Compute the brace value of a pulsar')
    parser.add_argument('-s', '--start_epoch', type=float,
                        help='Start epoch of observation (in MJD)', required=True)
    parser.add_argument('-e', '--end_epoch', type=float,
                         help='End epoch of observation (in MJD)', required=True)
    parser.add_argument('-p', '--parfile', type=str,
                        help='Parameter file', required=True)
    parser.add_argument("-t", "--telescope", type=str, help='Telescope name', default='MeerKAT')

    parser.add_argument("-f", "--fine", help = "fine search each conjunciton", action="store_true")

    parser.add_argument("-a", "--abs_ha", help = "HA to constrain the output", type=float) 

    parser.add_argument("-o",help='List of orbital phases that needs time')
    args = parser.parse_args()
    return args

def get_conjunction(par_file, start_epoch, end_epoch):
    pb_days = par_file.get_parameter('PB').value
    t0 = par_file.get_parameter('T0').value
    t0 = t0 + np.floor((start_epoch - t0)/pb_days)*pb_days
    hiT = start_epoch
    oldnu = calc_true_anom(hiT, par_file)
    it=0
    step = 1e-2
    while (1):
        hiT += step*pb_days
        newnu = calc_true_anom(hiT, par_file)
        if (newnu > np.pi/2 and oldnu < np.pi/2): break
        else: oldnu = newnu
        it = it + 1 
        if hiT > start_epoch + pb_days:
            print("No conjunctions found, trying the full range")
            loT = start_epoch
            hiT = loT + pb_days
            break
        # print("{:.5f} {:.5f} {:.5f}".format(hiT, newnu, oldnu))
        loT = hiT - step*pb_days

    #print("Searching between ", loT, " and ", hiT)

    conjT = bisect(solve_T_at_ntrue_anom, par_file, loT, hiT, TOL=1.0e-15, MAXIT=1000)
    if (conjT == 0.0):
        print("No conjunctions found")
        sys.exit(0)
    return conjT

def main():
    args = parse_args()
    par_file = ParFile(args.parfile)
    pb_days = par_file.get_parameter('PB').value
    start_epoch = args.start_epoch


    obs_location = EarthLocation.of_site("MeerKAT")
    observer = Observer(location=obs_location, name="MeerKAT", timezone="UTC")

    conjT = get_conjunction(par_file, args.start_epoch, args.start_epoch + pb_days)


    conjunctions = []
    if(args.fine):
        start_epochs = np.arange(conjT -0.1*pb_days, args.end_epoch, step=pb_days)
        for start_epoch in start_epochs:
            conjT = get_conjunction(par_file, start_epoch, start_epoch + pb_days)
            conjunctions.append(conjT)
    else:
        conjunctions = np.arange( conjT, args.end_epoch, step=pb_days)


    conjunctions = np.array(conjunctions)

    import pandas as pd
    df = None
    sky_coord = SkyCoord(par_file.get_parameter('RAJ').value, par_file.get_parameter('DECJ').value, unit=(u.hourangle, u.deg), frame='icrs')		

    for conjunction in conjunctions:
        time = Time(conjunction,format='mjd', scale='tdb')
        altaz_frame = AltAz(obstime=time, location=obs_location)
        source_altaz = sky_coord.transform_to(altaz_frame)
        prev_date = Time(time.mjd -1 , format='mjd')
        E1 = source_altaz.alt.deg    #
        Az = source_altaz.az.deg
        za = 90.0 - E1

        lst = observer.local_sidereal_time(time) 

        ha = (lst - sky_coord.ra).to(u.deg).value  # hour angle in degrees

        # Adjust HA to be within ±180 degrees if desired
        if ha > 180:
            ha -= 360
        elif ha < -180:
            ha += 360
        iso_time = str(time.iso).replace(" ","-")

        target = FixedTarget(coord=sky_coord, name=par_file.get_parameter('PSRJ').value)

        
        #pa = parallactic_angle(float(mkt.lat), float(ha), float(dec))
        try:
            rt = observer.target_rise_time(prev_date.iso, target, which="next", horizon=15*u.deg ).iso
            rt = str(rt).replace(" ","-").replace("/","-")

        except Exception as e:
            print (e)
            rt = 0
        try:
            st = observer.target_set_time(time.iso, target, which="next", horizon=15*u.deg ).iso
            st = str(st).replace(" ","-").replace("/","-")
        except Exception as e:
            print (e)
            st = 0
        topo_mjd = bary2topo(sky_coord, time)
        topo_utc = Time(topo_mjd, format='mjd', scale='utc')


        if df is None:
            df = pd.DataFrame([{'PSRJ': par_file.get_parameter('PSRJ').value, 'CONJUNCTION_TIME_MJD_BARY': time.mjd, 'CONJUNCTION_TIME_UTC_BARY': iso_time, 'CONJUNCTION_TIME_MJD_TOPO': topo_mjd, 'CONJUNCTION_TIME_UTC_TOPO': topo_utc.isot, 'HA': ha/15.0, 'LST': lst, 'ZA': za, 'RISE': rt, 'SET': st}])
        else:
            df = pd.concat([df, pd.DataFrame([{'PSRJ': par_file.get_parameter('PSRJ').value, 'CONJUNCTION_TIME_MJD_BARY': time.mjd, 'CONJUNCTION_TIME_UTC_BARY': iso_time, 'CONJUNCTION_TIME_MJD_TOPO': topo_mjd, 'CONJUNCTION_TIME_UTC_TOPO': topo_utc.isot, 'HA': ha/15.0, 'LST': lst, 'ZA': za, 'RISE': rt, 'SET': st}])], ignore_index=True)

        #print(f"{par_file.get_parameter('PSRJ').value} {time.mjd} {iso_time} {topo_mjd} {topo_utc.isot} {ha/15.0} {lst} {za} {rt} {st}")

    #print dataframe 

    #shortlist only the ones with HA within the range
    if args.abs_ha:
        df = df[abs(df['HA']) < args.abs_ha]

    pd.set_option('display.max_rows', None)

    df.drop_duplicates(inplace=True)
    #drop psrj, CONJUNCTION_TIME_MJD_BARY, CONJUNCTION_TIME_MJD_TOPO
    #df.drop(columns=['PSRJ', 'CONJUNCTION_TIME_MJD_BARY', 'CONJUNCTION_TIME_MJD_TOPO', 'CONJUNCTION_TIME_UTC_BARY'], inplace=True)

    #do not display index
    df.reset_index(drop=True, inplace=True)

    print(df)




if __name__ == "__main__":
    main()
