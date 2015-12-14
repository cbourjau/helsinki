import numpy as np

def read_mca(filename):
    '''A simple function to read .mca files
        - filename    path/filename to the .mca files
        returns
        data and live time
    '''
    with open(filename,'r') as f:
        lines = f.readlines()

        index = 0
        for l in lines:
            if("LIVE_TIME" in l):
                live_time = float(l.split()[2])
            if("<<DATA>>" in l):
                break
            index += 1
        data = list()
        for i in range(index+1,len(lines)-1):
            data.append(float(lines[i]))
    return  np.array(data), live_time


def pow(x,n,i):
    return n*x**i
def powegdi(x,n,i):
    return n*np.log(x)*x**i
def powegdn(x,n,i):
        return x**i
def lin_fun(x,k,m):
    return x*k+m


def calibration_params(gain, calibration):#:, slopes, slope_er, intercept, intercept_er):
    from scipy.optimize import curve_fit
    if(gain in calibration.keys()):
        slope = calibration[gain][0][0]
        intercept = calibration[gain][0][1]
        slope_err = calibration[gain][1][0][0]
        intercept_err = calibration[gain][1][1][1]
        return (slope,intercept),(slope_err,intercept_err)
    if(gain<100):
        raise ValueError("Don't know how to interpolate or extrapolate below 100 gain")
    slopes = list()
    intercepts = list()
    slopes_err = list()
    intercepts_err = list()

    for k,v in calibration.iteritems():
        slopes.append(v[0][0])
        intercepts.append(v[0][1])
        slopes_err.append(np.sqrt(v[1][0][0]))
        intercepts_err.append(np.sqrt(v[1][1][1]))

    pop, pcov = curve_fit(pow, calibration.keys(), slopes, p0 = [1.,-1.0], sigma = slopes_err)
    slope = pow(gain,pop[0],pop[1])


    err_slope = np.sqrt(powegdi(gain,pop[0],pop[1])**2*pcov[0][0] +
                        powegdn(gain,pop[0],pop[1])**2*pcov[1][1] +
                        2*powegdn(gain,pop[0],pop[1])*powegdi(gain,pop[0],pop[1])*pcov[1][0])
    intercept = calibration[100][0][1]
    err_intercept = calibration[100][0][1]/2
    return (slope,intercept),(err_slope,err_intercept)


def chantoV(chan, param, er_param):
    err = np.sqrt(chan**2*er_param[0]**2+er_param[1]**2)
    return (lin_fun(chan,param[0],param[1]),err)

def chantoE(chan, param, er_param):
    C = 1e13/1.60217662 #convertsion factor from mV to #electrons
    err = np.sqrt(chan**2*er_param[0]**2+er_param[1]**2)
    return (lin_fun(chan,param[0],param[1])*C,err*C)
