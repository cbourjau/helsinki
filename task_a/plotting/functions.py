#!/usr/bin/python
import numpy, scipy
from scipy import stats

def line(x, m, b):
    return m*x+b

def fit_slope_min(volts, caps, caps_errors, nvalues):
    npoints = len(volts)
    best_params = [1e9,1e9]
    best_errors = None
    start = int(npoints/2.)
    
    for i in range(start,npoints-nvalues):
        params, errors = scipy.optimize.curve_fit(line, volts[i:i+nvalues], caps[i:i+nvalues],
                                                  sigma = caps_errors[i:i+nvalues])
        if numpy.fabs(params[0]) < numpy.fabs(best_params[0]):
            best_params = params
            best_errors = errors

    return numpy.array(best_params), numpy.sqrt(numpy.diag(best_errors))

def fit_slope_max(volts, caps, caps_errors, nvalues):
    npoints = len(volts)
    best_params = [-1e9,-1e9]
    best_errors = None
    
    for i in range(npoints-nvalues):
        params, errors = scipy.optimize.curve_fit(line, volts[i:i+nvalues], caps[i:i+nvalues],
                                                  sigma = caps_errors[i:i+nvalues])
        if params[0] > best_params[0]:
            best_params = params
            best_errors = errors

    return numpy.array(best_params), numpy.sqrt(numpy.diag(best_errors))

def find_transition_point(p1, std1, p2, std2):

    vfd = (p2[1]-p1[1])/(p1[0]-p2[0])
    print "vfd", vfd#, p1, p2
    vfd_error = numpy.fabs((std2[1]-std1[1])/(p1[0]-p2[0]) - \
        vfd*(std1[0]-std2[0])/(p1[0]-p2[0]))
    print "vfd_error", vfd_error#, std1, std2
    
    return numpy.array([numpy.array(vfd), numpy.array(vfd_error)])
