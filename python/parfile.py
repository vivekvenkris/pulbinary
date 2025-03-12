import re


import constants
from constants import DEGREE_TO_RADIAN, fit_flags, str_keys, float_keys, aliases, floatn_keys, TSUN
from parameter import parameter
from exceptions import BinaryParameterMissing
from utils import convert_to_float, convert_to_int
import numpy as np
from uncertainties import umath,ufloat
from astropy.coordinates import SkyCoord



class ParFile(object):

    def __init__(self, file_name):
        self._file_name = file_name
        self._parameters = {}
        self._parse_par()
        self._check_and_consolidate_other_parameters()
        self._compute_mass_function()

    @property
    def file_name(self):
        return self._file_name

    @property
    def parameters(self):
        return self._parameters

    def get_parameter(self, param_name):
        if param_name in self._parameters:
            return self._parameters[param_name]

        else:
            raise ValueError("No parameter:" + param_name)

    def get_value(self, param_name):
        return self.get_parameter(param_name).value
    def get_error(self, param_name):
        return self.get_parameter(param_name).error        

    def has_parameter(self, param_name):
        return param_name in self._parameters



    def add_parameter(self, name, value, fit, error, alias=None):
        if name in self._parameters:
            raise ValueError("Parameter already exists:" + name)

        else:
            if alias is None:
                alias = aliases[name] if name in aliases else None

            p = parameter(name, value, fit, error, alias)
            self._parameters[name] = p

    def add_or_update_parameter(self, name, value, error, fit=None, alias=None):
        if name in self._parameters:
            self._parameters[name].update(name, value, error, fit, alias)

        else:
            self.add_parameter(name, value, fit, error, alias)

    def _add_parameter_from_ufloat(self, name, ufloat, fit=0, alias=None):
        self.add_parameter(name, ufloat.nominal_value, fit, ufloat.std_dev, alias)


    def _check_and_consolidate_other_parameters(self):

        self._parameters['TEMPO_VERSION'] = 1 \
            if 'T2CMETHOD' in self._parameters and self._parameters['T2CMETHOD'] == 'TEMPO' \
            else 2

        if 'BETA' in self._parameters and 'LAMBDA' in self._parameters:
            self._parameters['ELAT'] = self._parameters['BETA']
            self._parameters['ELON'] = self._parameters['LAMBDA']

        if 'BINARY' not in self._parameters:
            raise BinaryParameterMissing(parameter='BINARY')

        for p in ['PB', 'A1']:
            if p not in self._parameters:
                raise BinaryParameterMissing(parameter=p)

        if 'PB' in self._parameters:
            n_orbit = 2 * np.pi / (self._parameters['PB'].ve  * constants.SECONDS_IN_A_DAY)
            self._add_parameter_from_ufloat("N_ORBIT", n_orbit)
        if 'PSR' in self._parameters:
            self._parameters['PSRJ'] = self._parameters['PSR']
            

        if 'ELL1' not in self._parameters['BINARY'].value:

            if 'E' in self._parameters:
                self._parameters['TEMPO_VERSION'] = 1
                self._parameters['ECC'] = self._parameters['E']

            elif 'ECC' in self._parameters:
                self._parameters['TEMPO_VERSION'] = 2

            else:
                raise BinaryParameterMissing(parameter="E/ECC")

            for p in ['T0', 'OM']:
                if p not in self._parameters:
                    raise BinaryParameterMissing(parameter=p)
                
        else:
            for p in ['EPS1', 'EPS2', 'TASC']:
                if p not in self._parameters:
                    raise BinaryParameterMissing(parameter=p)
            ecc = umath.sqrt(self._parameters['EPS1'].ve ** 2 + self._parameters['EPS2'].ve ** 2)
            self._parameters['ECC'] = parameter('ECC', ecc.nominal_value, 0, ecc.std_dev,
                                               aliases['ECC'] if 'ECC' in aliases else None)
            om = umath.atan2(self._parameters['EPS1'].ve, self._parameters['EPS2'].ve)
            self._parameters['OM'] = parameter('OM', ecc.nominal_value, 0, ecc.std_dev,
                                              aliases['OM'] if 'OM' in aliases else None)
            self._parameters['T0'] = self._parameters['TASC']

        for p in ['PBDOT', 'XDOT', 'EDOT']:
            if p in self._parameters:
                if self._parameters['TEMPO_VERSION'] == 1 and self._parameters[p] is not None:
                    self._parameters[p].scale(1e-12)
        
        if self._parameters['TEMPO_VERSION'] == 1:
            if 'VARSIGMA' in self._parameters:
                self._parameters['SITG'] = self._parameters['VARSIGMA']

        if 'H3' in self._parameters and 'STIG' in self._parameters:
            m2, sini = self._shapiro_converter(direction="h3stig->m2sini")
            self._add_parameter_from_ufloat('M2', m2)
            self._add_parameter_from_ufloat('SINI', sini)
            self._add_parameter_from_ufloat('INC', umath.asin(sini)* 180.0/np.pi)

            
        elif 'M2' in self._parameters and 'SINI' in self._parameters:
            h3,stig = self._shapiro_converter(direction="m2sini->h3stig")
            self._add_parameter_from_ufloat('H3', h3)
            self._add_parameter_from_ufloat('STIG', stig)
            self._add_parameter_from_ufloat(
                'INC', umath.asin(self._parameters['SINI'].ve) * 180.0/np.pi)


        if 'PMRA' in self._parameters and 'PMDEC' in self._parameters:
            pm_total = umath.sqrt(self._parameters['PMRA'].ve ** 2 + self._parameters['PMDEC'].ve ** 2)
            self._add_parameter_from_ufloat('PM_TOTAL', pm_total)

        if 'P0' in self._parameters:
            f0 = 1/self._parameters['P0'].ve
            self._add_parameter_from_ufloat('F0', f0)

        elif 'F0' in self._parameters:
            p0 = 1/self._parameters['F0'].ve
            self._add_parameter_from_ufloat('P0', p0)

        if 'OMDOT' in self._parameters:
            m_from_omdot = ParFile.get_M_from_omdot(
                self._parameters['N_ORBIT'].ve, self._parameters['ECC'].ve, self._parameters['OMDOT'].ve)
            self._add_parameter_from_ufloat('M_FROM_OMDOT', m_from_omdot)


    @staticmethod
    def get_M_from_omdot(n_orbit, ecc, omdot):
        rest_omdot = 3.0 * \
                np.abs(constants.TSUN)**(2/3.0) * \
                n_orbit**(5/3.0) / (1-ecc*ecc)
        omdot_rad_per_sec = omdot * \
            constants.DEGREE_TO_RADIAN / constants.SECONDS_IN_A_YEAR
        total_mass = np.abs(omdot_rad_per_sec / rest_omdot) ** 1.5
        return total_mass

    @staticmethod
    def get_omdot_from_M(n_orbit, ecc, m):
        return 3.0 * n_orbit **(5/3.0) * np.abs(constants.TSUN * m)**(2/3.0) / (1-ecc*ecc) * constants.RADIAN_TO_DEGREE  * constants.SECONDS_IN_A_YEAR


    @staticmethod
    def get_from_mf(mf, m1, m2, sini):
        if m1 == "?":
            return np.sqrt(((m2 * sini) ** 3) / mf) - m2
        if m2 == "?":
            coeff = [sini ** 3, -1 * mf, -2 * mf * m1, -1 * mf * m1 * m1]
            roots = np.roots(coeff)
            roots = roots[np.isreal(roots)]
            return np.real(roots[0])
        
        mtot = m1 +m2
        if mf == "?":
            return ((m2 * sini) ** 3) / mtot ** 2
        if sini == "?":
            return (mf * mtot * mtot / (m2 ** 3)) ** 1/3.

    




    @staticmethod
    def get_h3_Stig(m2, sini):
        """
        Returns the orthometric parameters of the Shapiro delay
        """
        inc = umath.asin(sini)
        stig = sini / (1 + umath.cos(inc)) 
        h3 = m2 * TSUN * stig ** 3 
        return h3, stig

    @staticmethod
    def get_m2_sini(h3, stig):
        """
        Returns default Shapiro delay parameters from the orthometric parameters
        """
        # because stig = tan(i/2) so sini = sin(2*arctan(stig)) = 2stig/1+stig^2
        sini = 2 * stig / (1 + stig*stig)
        m2 = h3 / (stig ** 3 * TSUN)
        return m2, sini

    def _shapiro_converter(self, direction):
        if direction == "h3stig->m2sini":
            if 'H3' not in self._parameters or 'STIG' not in self._parameters:
                return None, None

            else:
                return ParFile.get_m2_sini(self._parameters['H3'].ve, self._parameters['STIG'].ve)

        elif direction == "m2sini->h3stig":
            if 'M2' not in self._parameters or 'SINI' not in self._parameters:
                return None, None

            else:
                return ParFile.get_h3_Stig(self._parameters['M2'].ve, self._parameters['SINI'].ve)


    def _compute_mass_function(self):
        if 'PB' not in self._parameters or 'A1' not in self._parameters:
            return

        else:
            pb = self._parameters['PB'].ve * constants.SECONDS_IN_A_DAY
            a1 = self._parameters['A1'].ve

            mass_func = 4 * np.pi * np.pi * a1 ** 3 / (pb ** 2 * constants.TSUN)

            p = parameter("MASS_FUNC", mass_func.n, 0, mass_func.s, None)
            self._parameters["MASS_FUNC"] = p

    def _parse_par(self):

        lines = open(self.file_name).readlines()
        for line in lines:

            param = value = error = fit = alias = None
            # skip if the line is empty or a comment
            if not line or '#' == line[0] or 'JUMP' in line or 'TNEF' in line or 'TNEQ' in line or line.strip() == "":
                continue

            # Convert any 'D-' or 'D+' to 'E-' or 'E+', strip leading or trailing space, and split
            chunks = line.replace("D-", "E-").replace("D+", "E+").strip().split()
            num_chunks = len(chunks)
            key = chunks[0]

            param = key

            if key in constants.ignore_keys:
                continue

            # Regex checks for non-digit chars, followed by digit chars
            m1 = re.search(r'(\D+)(\d+)$', param)

            # This one looks for the DMX[RF][12]_* params
            m2 = re.search(r'(\D+\d+_)(\d+)$', key)

            m = m1 if m1 is not None else m2

            if key == "JUMP":
                if num_chunks > 3:

                    param = param + '_%s_%s' % (chunks[1], chunks[2])

                    if chunks[3] not in fit_flags:
                        value = convert_to_float(chunks[3])

                        if num_chunks == 5:  # Eg: JUMP -be MK -0.000306243 0

                            if chunks[4] not in fit_flags:  # Eg: JUMP -be MK -0.000306243 0.01
                                error = convert_to_float(chunks[4])

                            elif num_chunks == 6:  # Eg: JUMP -be MK -0.000306243 0 0.01
                                fit = convert_to_int(chunks[4])
                                error = convert_to_float(chunks[5])

            elif key in str_keys:
                value = chunks[1]


            elif key in float_keys:
                if param == 'SINI' and chunks[1] == 'KIN':
                    print("SINI is KIN")
                    continue
                value = convert_to_float(chunks[1])


            elif m is not None and m.group(1) in floatn_keys:
                value = convert_to_float(chunks[1])

            if num_chunks == 3 and chunks[2] not in fit_flags:  # Eg: F0 1.2 0.000002
                error = convert_to_float(chunks[2])

            elif num_chunks == 4:  # Eg: F0 1.2  1 0.000002
                fit = convert_to_int(chunks[2])
                error = convert_to_float(chunks[3])

            p = parameter(param, value, fit, error,
                          aliases[param] if param in aliases else None)
            self._parameters[param] = p

            if param == 'KIN':
                self._parameters['SINI'] = parameter('SINI', umath.sin(value * DEGREE_TO_RADIAN), fit, umath.sin(error * DEGREE_TO_RADIAN) if error is not None else None,
                          aliases[param] if param in aliases else None)


