import numpy as np
import scipy
from scipy.optimize import curve_fit

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
def exp_fit(x,a,k): return a*np.exp(k*x)
def gauss_fit(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))
def double_gauss_fit(x, a1, x1, sigma1, a2, x2, sigma2):
    return gauss_fit(x,a1,x1,sigma1) + gauss_fit(x,a2,x2,sigma2)

####################################################################
# Derivatives of the interpolation fits
####################################################################
def dConst_fit(x,k, sk): return 
def dLinear_fit(x,m,b,sm,sb): return 
def dPower_fit(x,m,k,sm,sk): return 
def dExp_fit(x,a,k): return a*k*np.exp(k*x)
def dGauss_fit(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2)) * (-2*(x-x0)/2*sigma**2)
def dDouble_gauss_fit(x, a1, x1, sigma1, a2, x2, sigma2):
    return dGauss_fit(x,a1,x1,sigma1) + dGauss_fit(x,a2,x2,sigma2)

####################################################################
# Create the interpolation for the gains
####################################################################
def fitGainComponents(calibration={}):
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
    sigma_slope = np.sqrt(np.diagonal(slope_cov)
    dSlope_d = 

    return charges

####################################################################
# Fit the noise data to an exponential rate vs channel number
####################################################################
def fitNoise(noise_data={}):
    noise_fits = {}
    noise_cuts = {}

    for gain in noise_data.keys():
        # find the bin with actual data!
        max_bin = 0
        max_rate = 0
        first_zero_bin = -1
        for bin in range(len(noise_data[gain][0])):
            if noise_data[gain][1][bin] > max_rate: 
                max_bin = bin
                max_rate = noise_data[gain][1][bin]
            if noise_data[gain][1][bin] == 0 and not max_bin==0:
                first_zero_bin = bin
                break

        params, covariance = curve_fit(linear_fit,
                                       noise_data[gain][0][max_bin:first_zero_bin],
                                       np.log(noise_data[gain][1][max_bin:first_zero_bin]),
                                       #p0=[-2e-1, 3000],
                                       )
        noise_fits[gain] = [params, covariance]

        # Find the bin where the noise decreases to less than 5% of the total
        total_rate = np.sum(noise_data[gain][1])
        noise_cuts[gain] = noise_data[gain][0][max_bin:][noise_data[gain][1][max_bin:] < 0.01*total_rate][0]

    return noise_fits, noise_cuts
    
####################################################################
# Fit a double gaussian to each
####################################################################
def fitDoubleGaussians(data={}):
    fitparams = {}
    
    voltages = data.keys()
    voltages.sort()
    for voltage in voltages:

        # some simple peakfinding to find the first and second peaks
        peak2_value = np.max(data[voltage][1])
        peak2_charge = [data[voltage][0][x] for x in range(len(data[voltage][0])) if data[voltage][1][x]==peak2_value][0]
        peak2_bin = [x for x in range(len(data[voltage][0])) if data[voltage][1][x]==peak2_value][0]

        # some simple peakfinding to find the first and second peaks
        peak1_value = np.max(data[voltage][1][:(peak2_bin*2)/3])
        peak1_charge = [data[voltage][0][x] for x in range(len(data[voltage][0][:(peak2_bin*2)/3])) if data[voltage][1][x]==peak1_value][0]
        peak1_bin = [x for x in range(len(data[voltage][0])) if data[voltage][1][x]==peak1_value][0]

        # Preliminary fits to get better estimates of parameters...
        params1, covariance = curve_fit(gauss_fit,
                                       data[voltage][0][peak1_bin-30:peak1_bin+30],
                                       data[voltage][1][peak1_bin-30:peak1_bin+30],
                                       p0 = [peak1_value, peak1_charge, 0.1*peak1_charge])

        params2, covariance = curve_fit(gauss_fit,
                                       data[voltage][0][peak2_bin-30:peak2_bin+30],
                                       data[voltage][1][peak2_bin-30:peak2_bin+30],
                                       p0 = [peak2_value, peak2_charge, 0.1*peak2_charge])

        # Now fit both gaussians simultaneously
        params, covariance = curve_fit(double_gauss_fit,
                                       data[voltage][0],
                                       data[voltage][1],
                                       p0 = [params1[0], params1[1], params1[2],
                                             params2[0], params2[1], params2[2]])
        
        fitparams[voltage] = [params, covariance]

    return fitparams
