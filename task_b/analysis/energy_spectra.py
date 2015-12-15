#!/usr/bin/python
import os, sys, glob, pickle, numpy, scipy, matplotlib
from scipy.optimize import curve_fit
from matplotlib import pyplot
from utils import *

####################################################################
#
# Calibration of the gains
# Read the data for the calibration and extrapolate/interpolate
#
####################################################################
calibration_file = os.path.abspath("../data/calibration.pkl")
if not os.path.exists(calibration_file):
    print "No calibration file found. Run the calibration.py script in " \
        "this folder before continuing."
    sys.exit()

calib = pickle.load(open(calibration_file))

slope_params, slope_covariance, intercept_params, intercept_covariance = fitGainComponents(calib)

print numpy.sqrt(numpy.diagonal(slope_covariance))
print numpy.sqrt(numpy.diagonal(intercept_covariance))

print slope_params
print intercept_params

####################################################################
#
# Noise analysis
# Need this to subtract out the noise spectrum from our data later
#
####################################################################
noise_scan = {}
noise_fits = {}
for filepath in glob.glob(os.path.abspath("../data/noise*.mca")):
    fname = os.path.basename(filepath)
    print fname
    voltage = float(fname.split("_")[0].replace("noiseHV",""))
    gain = int(fname.split("_")[-1][4:-4])
    data, livetime = read_mca(filepath)

    # smooth over 3 bins
    #origdata = data[:]
    #for i in range(len(data[4:-4])):
    #    data[i] = numpy.sum(origdata[i-4:i+4])/9.0
    #print data

    charges = apply_calibration(data, slope_params, intercept_params, gain)
    # Save the noise as channel number, error, rate, error
    noise_scan[gain] = [numpy.cumsum(numpy.ones_like(data)), data/livetime, numpy.sqrt(data)/livetime]
    
noise_fits, cuts = fitNoise(noise_scan)

# Plot the spectra
pyplot.figure()
#pyplot.errorbar(noise_scan[50][0], noise_scan[50][2],
#                0,0, #noise_scan[50][3], 0,#noise_scan[50][3],                
#                label="Gain = %i" % 50, color='r')
#pyplot.errorbar(noise_scan[200][0], noise_scan[200][2],
#                0,0,#noise_scan[200][3], 0,#noise_scan[200][3],                
#                label="Gain = %i" % 200, color='g')
#pyplot.errorbar(noise_scan[500][0], noise_scan[500][2],
#                0,0,#noise_scan[500][3], 0,#noise_scan[500][3],                
#                label="Gain = %i" % 500, color='b')


pyplot.errorbar(noise_scan[50][0], noise_scan[50][2]/numpy.sum(noise_scan[50][2]), 
                0,0, #noise_scan[50][3], 0,#noise_scan[50][3],                
                label="Gain = %i" % 50, color='r')
pyplot.errorbar(noise_scan[200][0], noise_scan[200][2]/numpy.sum(noise_scan[200][2]), 
                0,0,#noise_scan[200][3], 0,#noise_scan[200][3],                
                label="Gain = %i" % 200, color='g')
pyplot.errorbar(noise_scan[500][0], noise_scan[500][2]/numpy.sum(noise_scan[500][2]), 
                0,0,#noise_scan[500][3], 0,#noise_scan[500][3],                
                label="Gain = %i" % 500, color='b')

pyplot.legend()
pyplot.grid()
#pyplot.yscale('log')
pyplot.xlim(min(noise_scan[gain][0]), max(noise_scan[gain][0]))
#pyplot.ylim(min(noise_scan[gain][2]), max(noise_scan[gain][2]))
#pyplot.xlim(0, 75)
pyplot.xlabel("Charge (pC)")
pyplot.ylabel("Rate (Hz)")
pyplot.savefig("noises.pdf")
    


####################################################################
#
# Iron data with the voltage scan
#
####################################################################
voltage_scan = {}
voltage_params = {}
files = glob.glob(os.path.abspath("../data/voltage_scan*.mca"))
files.sort()
f, axes = pyplot.subplots(4,3)#, sharex=True, sharey=True)
i=0
j=0

for filepath in files: 
    fname = os.path.basename(filepath)
    voltage = float(fname.split("_")[2].replace("HV",""))
    gain = float(fname.split("_")[-1][4:-4])

    # Get the data
    data, livetime = read_mca(filepath)
    data = numpy.array(data[45:])

    # Apply the calibration using the gain
    # Get this slope and intercept
    charges = apply_calibration(data, slope_params, intercept_params, gain)
    voltage_scan[voltage] = [charges, data/livetime, numpy.sqrt(data)/livetime]
        
    fit = fitDoubleGaussians({voltage:[charges*1e12, data/livetime, numpy.sqrt(data)/livetime]})[voltage]
    voltage_params[voltage] = fit
    axes[i,j].scatter(charges*1e12, data/livetime)
    axes[i,j].set_xlim(numpy.min(charges*1e12), numpy.max(charges*1e12))
    axes[i,j].set_ylim(numpy.min(data/livetime), numpy.max(data/livetime))

    a1 = fit[0][0]
    x1 = fit[0][1]
    s1 = fit[0][2]
    a2 = fit[0][3]
    x2 = fit[0][4]
    s2 = fit[0][5]
    
    fitvalues = double_gauss_fit(charges*1e12, a1, x1, s1, a2, x2, s2)
    axes[i,j].plot(charges*1e12, fitvalues, 'rx')    
    
    i += 1
    if i == 4:
        j += 1
        i = 0

pyplot.tight_layout()
pyplot.savefig("voltage_scan.pdf")

# Now we have values. Hurrah!
# lets plot the peak positions vs the voltage just for fun
voltages = voltage_params.keys()
voltages.sort()
peak1 = [voltage_params[v][0][1] for v in voltages]
peak2 = [voltage_params[v][0][4] for v in voltages]

pyplot.figure()
pyplot.scatter(voltages, peak1, color='r')
pyplot.scatter(voltages, peak2, color='b')
pyplot.yscale("log")
pyplot.grid()
pyplot.xlabel("High Voltage (kV)")
pyplot.ylabel("Peak Position (pC)")
pyplot.savefig("peak_positions.pdf")

# And do the resolution (ie, sigma) vs voltage
sigma1 = [voltage_params[v][0][2] for v in voltages]
sigma2 = [voltage_params[v][0][5] for v in voltages]

pyplot.figure()
pyplot.scatter(voltages, sigma1, color='r')
pyplot.scatter(voltages, sigma2, color='b')
pyplot.yscale("log")
pyplot.grid()
pyplot.xlabel("High Voltage (kV)")
pyplot.ylabel("Charge Resolution")
pyplot.savefig("peak_resolutions.pdf")

sys.exit()


####################################################################
#
# Cesium analysis
# These values are mV/channel and are for each gain in the form
# of a {gain: [slope, intercept]} from a plot of the mV vs channel
#
####################################################################
ce_data, ce_livetime = read_mca(os.path.abspath("../data/cesiumHV0.96_gain50.mca"))
ce_gain = 50

ce_slope = slope_slope * ce_gain + slope_intercept
ce_intercept = intercept_slope * ce_gain + intercept_intercept

# Get the charge for this data
ce_charges = numpy.cumsum(numpy.ones_like(ce_data)) # Channel numbers
ce_charges = ce_slope * ce_charges + ce_intercept # Now in mV
ce_charges *= mVolts * amplifier_capacitance # Now in Coulombs

print ce_charges
pyplot.figure()
pyplot.scatter(ce_charges/pC, ce_data/(ce_livetime))
pyplot.grid()
pyplot.xlim(min(ce_charges/pC), max(ce_charges/pC))
pyplot.xlabel("Charge (pC)")
pyplot.ylabel("Rate (Hz)")
pyplot.savefig("cesium.pdf")
