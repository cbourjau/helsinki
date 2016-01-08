import numpy as np
import scipy
from scipy.optimize import curve_fit
from scipy.odr import ODR, Model, Data, RealData

mVolts = 1e-3
pC = 1e-12
amplifier_capacitance = 1000.0e-12 # Farads

####################################################################
# Read the data from the .mca file and return less useless stuff
####################################################################
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

####################################################################
# fits for interpolation
####################################################################
def const_fit(x,k): return k
def linear_fit(x,m,b): return m*x+b
def power_fit(x,m,k): return m*x**k
def exp_fit(x,a,k): return a*np.exp(-k*x)
def exp2_fit(x,a,k,c): return a*np.exp(-k*x)+c
def gauss_fit2(x,a,x0,s): return gauss_fit([a,x0,s],x)

def gauss_fit(params, x):
    if len(params) < 3: return 1e10
    return params[0]*np.exp(-(x-params[1])**2/(2*params[2]**2))
def double_gauss_fit(params, x):
    if len(params) < 6: return 1e10
    return gauss_fit(params[:3],x) + gauss_fit(params[3:],x)

def powegdi(x,n,i): return n*np.log(x)*x**i
def powegdn(x,n,i): return x**i


def chantoV(chan, param, er_param):
    C = amplifier_capacitance
    err = np.sqrt(chan**2*er_param[0]**2+er_param[1]**2)
    return (linear_fit(chan,param[0],param[1])*C,err*C)

def chantoE(chan, param, er_param):
    C = amplifier_capacitance/1.60217662e-19 #conversion factor from mV to #electrons
    err = np.sqrt(chan**2*er_param[0]**2+er_param[1]**2)
    return (linear_fit(chan,param[0],param[1])*C,err*C)



####################################################################
# Create the interpolation for the gains
####################################################################
def calibration_params(gain, calibration):#:, slopes, slope_er, intercept, intercept_er):

    if(gain in calibration.keys()):
        slope = calibration[gain][0][0]
        intercept = calibration[gain][0][1]
        slope_err = calibration[gain][1][0][0]
        intercept_err = calibration[gain][1][1][1]
        return (slope,intercept),(slope_err,intercept_err)

    if(gain<20):
        raise ValueError("Don't know how to interpolate or extrapolate below 20 gain")
    slopes = list()
    intercepts = list()
    slopes_err = list()
    intercepts_err = list()

    for k,v in calibration.iteritems():
        slopes.append(v[0][0])
        intercepts.append(v[0][1])
        slopes_err.append(np.sqrt(v[1][0][0]))
        intercepts_err.append(np.sqrt(v[1][1][1]))


    pop, pcov = curve_fit(power_fit, calibration.keys(), slopes, p0 = [1.,-1.0], sigma = slopes_err)
    slope = power_fit(gain,pop[0],pop[1])

    err_slope = np.sqrt(powegdi(gain,pop[0],pop[1])**2*pcov[0][0] +
                        powegdn(gain,pop[0],pop[1])**2*pcov[1][1] +
                        2*powegdn(gain,pop[0],pop[1])*powegdi(gain,pop[0],pop[1])*pcov[1][0])
    intercept = calibration[100][0][1]
    err_intercept = calibration[100][0][1]/2.0#TODO: is this the right factor?
    return (slope,intercept),(err_slope,err_intercept)

    gains = []
    slopes = []
    intercepts = []
    for gain in calibration.keys():

        gains.append(float(gain.replace("gain","")))
        slopes.append(calibration[gain][0][0])
        intercepts.append(calibration[gain][0][1])
        cov = calibration[gain][1]

    popt_slopes, pcov_slopes = curve_fit(power_fit, gains, slopes, sigma=np.sqrt(np.diagonal(cov))[0] )
    popt_intercepts, pcov_intercepts = curve_fit(const_fit, gains, intercepts, sigma=np.sqrt(np.diagonal(cov))[1] )
    #popt_intercepts = np.mean(intercepts)

    return popt_slopes, pcov_slopes, popt_intercepts, pcov_intercepts

####################################################################
# Apply the calibration using fits to the slope and intercept
# values from Samuel's fit
####################################################################
def apply_calibration(data,
                      slope_params, slope_cov,
                      intercept_params, intercept_cov,
                      gain, capacitance=1000.0e-12):
    slope = power_fit(gain, slope_params[0], slope_params[1])
    intercept = const_fit(gain, intercept_params[0])

    # Get the charge for this data
    charges = np.cumsum(np.ones_like(data)) # Channel numbers
    charges = slope * charges + intercept # Now in mV
    charges *= mVolts * capacitance # Now in Coulombs

    # Charge errors?
    # slope error first
    sigma_slope = np.sqrt(np.diagonal(slope_cov))

    return charges

####################################################################
# Fit the noise data to an exponential rate vs channel number
####################################################################
def fitNoise(noise_data={}):
    noise_fits = {}
    noise_cuts = {}

    # find the bin with actual data!
    max_bin = 0
    max_rate = 0
    first_zero_bin = -1
    for bin in range(len(noise_data[0])):
        if noise_data[1][bin] > max_rate:
            max_bin = bin
            max_rate = noise_data[1][bin]
        if noise_data[1][bin] == 0 and not max_bin==0:
            first_zero_bin = bin
            break

    #params, covariance = curve_fit(linear_fit,
    #                               noise_data[0][max_bin:first_zero_bin],
    #                               np.log(noise_data[1][max_bin:first_zero_bin]),
    #                               p0=[-2e-1, -0.5, 1e-2],
    #                               )

    params, covariance = curve_fit(exp2_fit,
                                   noise_data[0][max_bin:first_zero_bin],
                                   noise_data[1][max_bin:first_zero_bin],
                                   p0=[0.25, 0.075, 0.005],
                                   )

    noise_fits = [params, covariance]
    print noise_fits[0]

    # Find the bin where the noise decreases to less than 5% of the total
    total_rate = np.sum(noise_data[1])
    noise_cuts = noise_data[0][max_bin:][noise_data[1][max_bin:] < 0.01*total_rate][0]

    return noise_fits, noise_cuts

####################################################################
# Fit a double gaussian to each
####################################################################
def fitDoubleGaussians(charges, charge_err, data, data_err):

    # some simple peakfinding to find the first and second peaks
    peak1_value = np.max(data)
    peak1_bin = [x for x in range(len(data)) if data[x]==peak1_value][0]
    peak1_charge = charges[peak1_bin]

    # Preliminary fits to get better estimates of parameters...
    first_peak_data = RealData(charges[peak1_bin-30:peak1_bin+30],
                               data[peak1_bin-30:peak1_bin+30],
                               charge_err[peak1_bin-30:peak1_bin+30],
                               data_err[peak1_bin-30:peak1_bin+30],)
    first_peak_model = Model(gauss_fit)
    first_peak_odr = ODR(first_peak_data, first_peak_model,
                         [peak1_value, peak1_charge, 0.1*peak1_charge])
    first_peak_odr.set_job(fit_type=2)
    first_peak_output = first_peak_odr.run()
    first_peak_params = first_peak_output.beta

    #second_peak_params, covariance = curve_fit(gauss_fit2,
    #                                           data[peak1_bin-30:peak1_bin+30],
    #                                           data[peak1_bin-30:peak1_bin+30],
    #                                           p0 = [peak1_value, peak1_charge, 0.1*peak1_charge])


    # subtract the largest peak so we can search for the other one
    updated_data = data-gauss_fit(first_peak_params, charges)
    #updated_data = data[:int(len(data)*2.0/3.0)]

    peak2_value = np.max(updated_data)
    peak2_bin = [x for x in range(len(updated_data)) if updated_data[x]==peak2_value][0]
    peak2_charge = charges[peak2_bin]

    #first_peak_params, covariance = curve_fit(gauss_fit2,
    #                                          updated_data[peak2_bin-30:peak2_bin+30],
    #                                          updated_data[peak2_bin-30:peak2_bin+30],
    #                                          p0 = [peak2_value, peak2_charge, 0.1*peak2_charge])

    # and the second peak...
    second_peak_data = RealData(charges[peak2_bin-30:peak2_bin+30],
                               data[peak2_bin-30:peak2_bin+30],
                               charge_err[peak2_bin-30:peak2_bin+30],
                               data_err[peak2_bin-30:peak2_bin+30],)
    second_peak_model = Model(gauss_fit)
    second_peak_odr = ODR(second_peak_data, second_peak_model,
                         [peak2_value, peak2_charge, first_peak_params[2]])
    second_peak_odr.set_job(fit_type=2)
    second_peak_output = second_peak_odr.run()
    second_peak_params = second_peak_output.beta

    # Now fit both gaussians simultaneously
    double_peak_data = RealData(charges, data,
                                charge_err, data_err)
    double_peak_model = Model(double_gauss_fit)
    double_peak_odr = ODR(double_peak_data, double_peak_model,
                          [first_peak_params[0], first_peak_params[1], first_peak_params[2],
                           second_peak_params[0], second_peak_params[1], second_peak_params[2]])
    double_peak_odr.set_job(fit_type=0)
    double_peak_output = double_peak_odr.run()
    double_peak_params = double_peak_output.beta
    double_peak_sigmas = double_peak_output.sd_beta

    #double_peak_params = [first_peak_params[0], first_peak_params[1], first_peak_params[2],
    #                       second_peak_params[0], second_peak_params[1], second_peak_params[2]]

    return double_peak_params, double_peak_sigmas
