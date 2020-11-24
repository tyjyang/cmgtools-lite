import os
import subprocess
import fnmatch
import math
import numpy as np
import itertools

def relative_round(value, relative_digits):
    """Rounds to a given relative precision"""

    if value == 0 or isinstance(value, str) or np.isnan(value) or np.isinf(value):
        return value

    if isinstance(value, tuple):
        return (relative_round(x, relative_digits) for x in value)

    precision = math.ceil(math.log10(abs(value)))

    absolute_digits = - precision + relative_digits

    return round(value, int(absolute_digits))


def getNumberPrecision(value):

    if value == 0 or isinstance(value, str) or np.isnan(value) or np.isinf(value):
        return value
    if isinstance(value, tuple):
        return (getNumberPrecision(x) for x in value)

    return math.ceil(math.log10(abs(value)))


def getValuePrecisionWrtUncertainty(value,uncertainty):
    return getNumberPrecision(value) - getNumberPrecision(uncertainty)


def roundValueAndUncertainty(cont, valKey="y", uncKey="dy", sigDigitsUnc=2):

    #cont[uncKey] = [relative_round(x,sigDigitsUnc) for x in cont[uncKey]]
    #cont[valKey] = [relative_round(x,int(sigDigitsUnc+getValuePrecisionWrtUncertainty(x,e))) for x,e in itertools.izip(cont[valKey],cont[uncKey])]

    # maybe faster looping on elements
    for i,(v,e) in enumerate(itertools.izip(cont[valKey],cont[uncKey])):
        cont[uncKey][i] = relative_round(e,sigDigitsUnc)
        cont[valKey][i] = relative_round(v,int(sigDigitsUnc+getValuePrecisionWrtUncertainty(v,e)))


if __name__ == "__main__":


    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options]')
    parser.add_option('--sig-digits',
                      dest='sigDigits', 
                      type='int', default='2',
                      help='Number of significant digits to round uncertainty')
    (options,args) = parser.parse_args()


    # print relative_round(1.3456,2)    
    # print relative_round(0.3456,2)
    # print relative_round(912.3456,2)
    # print relative_round(0.0456,2)

    # print getNumberPrecision(1.3456)    
    # print getNumberPrecision(0.3456)
    # print getNumberPrecision(912.3456)
    # print getNumberPrecision(0.0456)

    sigDigits = options.sigDigits # for uncertainty, central to be rounded accordingly
    pairs = [(137.65, 0.35), 
             (0.00835,0.00012), 
             (170.5,80.3), 
             (0.00456,0.00165),
             (0.01,0.35),
             (16789, 1234)]
    print "before rounding\t\t after rounding"
    for p in pairs:
        before = str(p)
        v,e = p
        e = relative_round(e,sigDigits)
        v = relative_round(v,int(sigDigits+getValuePrecisionWrtUncertainty(v,e)))
        newp = (v,e)
        after = str(newp)
        print "{b}\t\t {a}".format(b=before,a=after)
