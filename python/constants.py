import numpy as np
from uncertainties import ufloat
import astropy.units as u


float_keys = ["PEPOCH", "POSEPOCH", "DM", "START", "FINISH",
              "TRES", "TZRMJD", "TZRFRQ", "NITS",
              "XDOT", "E", "ECC", "EDOT", "T0", "PB", "PBDOT", "OM", "OMDOT",
              "EPS1", "EPS2", "EPS1DOT", "EPS2DOT", "TASC", "LAMBDA", "BETA",
              "RA_RAD", "DEC_RAD", "GAMMA", "SINI", "M2", "MTOT", "XPBDOT",
              "ELAT", "ELONG", "PMLAMBDA", "PMBETA", "PX", "PMRA", "PMDEC", "OM2DOT", "X2DOT"
              "PB_2", "A1_2", "E_2", "T0_2", "OM_2", "DMX", "SOLARN0", "MODE", "NPRNT", "DR", "DTHETA", "H3", "STIG", "KIN", "KOM"]

floatn_keys = ["F", "P", "FB", "A", "FD", "DMX_", "DMXEP_", "DMXR1_",
               "DMXR2_", "DMXF1_", "DMXF2_"]

str_keys = ["FILE", "PSR", "PSRJ", "RAJ", "DECJ", "EPHEM", "CLK", "BINARY",
            "TZRSITE", "UNITS", "TIMEEPH", "T2CMETHOD", "CORRECT_TROPOSPHERE", "PLANET_SHAPIRO", "DILATEFREQ",
            "TZRSITE"]

ignore_keys= ["CHI2R", "NTOA", "NITS"]
fit_flags = ['0', '1']

aliases = {
    "ECC": "E",
    "E":  "ECC"

}

# MJD of the J2000.0 epoch
J2000 = 51544.5
TSUN = 0.000004925490947
SECONDS_IN_A_DAY = 86400.0
DEGREE_TO_RADIAN = np.pi / 180.0
RADIAN_TO_DEGREE = 1/DEGREE_TO_RADIAN
DAYS_IN_A_YEAR = 365.2425
SECONDS_IN_A_YEAR = SECONDS_IN_A_DAY * DAYS_IN_A_YEAR
GRAVITATIONAL_CONSTANT = 6.67408e-11
GMSUN = 1.32712438e+20
MSUN = 1.989e+30
SPEED_OF_LIGHT = 299792458.0
MILLI_ARC_SECCONDS_TO_DEGREE = 1 / 60 / 60 / 1000
KPC_TO_M= 3.08567758e19

#Gravity collaboration 2019
SSB_THETA_0 = ufloat(236.9,4.2) * 1000 # m s^-1
SSB_R_0 = ufloat(8.178, 0.026) * KPC_TO_M  # m

