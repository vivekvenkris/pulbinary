
import constants
from uncertainties import ufloat
import utils


class parameter(object):

    def __init__(self, name, value, fit, error, alias):
        self._name = name
        self._value = value
        self._fit = fit
        self._error = error
        self._alias = alias
        self._ve = ufloat(self._value, self._error if self._error is not None else 0) if utils.is_number(self._value) else None


    @property
    def name(self):
        return self._name

    @property
    def value(self):
        return self._value

    @property
    def fit(self):
        return self._fit

    @property
    def error(self):
        return self._error

    @property
    def alias(self):
        return self._alias

    @property
    def ve(self):
        return ufloat(self._value, self._error if self._error is not None else 0) if utils.is_number(self._value) else None
    
    def update(self, name =None, value=None, fit=None, error=None, alias=None):
        if name is not None:
            self._name = name

        if value is not None:
            self._value = value
        
        if fit is not None:
            self._fit = fit
        
        if error is not None:
            self._error = error
        
        if alias is not None:
            self._alias = alias



    def __repr__(self):
        return self.__str__()

    def scale(self, scale_value):
        if self._value is not None:
            self._value *= scale_value
        if self._error is not None:
            self._error *= scale_value

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.name == other.name
        else:
            return False

    def __cmp__(self, other):
        if isinstance(other, self.__class__):
            return self.name in other.name
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return self.name

    def __str__(self):
        return "Name:" + str(self.name) + " Value: " + str(self.value) + " Fit:  " + str(self.fit) + " Error: " \
               + str(self.error) + " Alias: " + str(self.alias) + "ve:" + str(self.ve)


