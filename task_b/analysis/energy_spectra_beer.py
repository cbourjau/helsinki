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

####################################################################
#
# Noise analysis
# Need this to subtract out the noise spectrum from our data later
#
####################################################################
noise_scan = {}
noise_fits = {}
print glob.glob(os.path.abspath("../data/beer can/group 1/*"))
for filepath in glob.glob(os.path.abspath("../data/beer can/group 1/Beercan_background_1600V_G100.mca")):
    fname = os.path.basename(filepath)
    print fname
    voltage = 1600
    gain = 100
    data, livetime = read_mca(filepath)
    bins = numpy.cumsum(numpy.ones_like(data))

    #data = data[10:60]
    #bins = bins[10:60]

    newbins = []
    newdata = []
    for bin in numpy.arange(0,len(bins),5):
        newdata.append(numpy.sum(data[bin:bin+4]))
        newbins.append(bins[bin])

    bins = numpy.array(newbins)
    data = numpy.array(newdata)
    
    if gain < 20: continue

    # Save the noise as channel number, error, rate, error
    noise_scan = [bins, data/livetime, numpy.sqrt(data)/livetime]

    #fit it?
    noise_fits, noise_cuts = fitNoise(noise_scan)


# Plot the spectra
pyplot.figure()
pyplot.errorbar(noise_scan[0], noise_scan[2],
                0,0,
                label="Gain = %i" % 100, color='r')
pyplot.scatter(noise_scan[0], exp2_fit(noise_scan[0], noise_fits[0][0], noise_fits[0][1], noise_fits[0][2]),
               label="Exponential Fit", color='b')
#pyplot.ylim(0,0.2)
pyplot.legend()
pyplot.grid()
pyplot.xlim(min(noise_scan[0]), max(noise_scan[0]))
pyplot.xlabel("Charge (pC)")
pyplot.ylabel("Rate (Hz)")
pyplot.savefig("beercan_noise.pdf")
    


####################################################################
#
# Iron data with the voltage scan
#
####################################################################
voltage_scan = {}
voltage_params = {}
files = glob.glob(os.path.abspath("../data/beer can/group 1/Beercan_Fe*.mca"))
files.sort()
f, axes = pyplot.subplots(2,2)#, sharex=True, sharey=True)
i=0
j=0

for filepath in files: 
    fname = os.path.basename(filepath)
    voltage = float(fname.split("_")[2].replace("V",""))
    gain = float(fname.split("_")[-1].split(".")[0][1:])

    if gain < 20: continue

    # Get the data
    data, livetime = read_mca(filepath)
    #data = numpy.array(data[45:]) # Removing channels below 45 since that seems to effectively remove noise

    #print numpy.sum(data/livetime)
    # Apply the calibration using the gain
    # Get this slope and intercept
    fitparams, fiterrors = calibration_params(gain, calib)
    
    channels = numpy.cumsum(numpy.ones_like(data))
    mvolts, charge_err = chantoV(channels, fitparams, fiterrors)

    voltage_scan[voltage] = [mvolts, charge_err, data/livetime, numpy.sqrt(data)/livetime]
    #fit, sigma = fitDoubleGaussians(mvolts, charge_err, data/livetime, numpy.sqrt(data)/livetime)

    #voltage_params[voltage] = [fit, sigma]
    axes[i,j].errorbar(x=mvolts, y=data/livetime, xerr=charge_err, yerr=numpy.sqrt(data)/livetime, label="Data")
    axes[i,j].set_xlim(numpy.min(mvolts), numpy.max(mvolts))
    axes[i,j].set_ylim(numpy.min(data/livetime), numpy.max(data/livetime))

    #fitvalues = double_gauss_fit(fit, mvolts)
    #axes[i,j].plot(mvolts, fitvalues, 'rx', label="Fit")
    
    i += 1
    if i == 2:
        j += 1
        i = 0

pyplot.tight_layout()
pyplot.xlabel("Charge (e)")
pyplot.ylabel("Events")
pyplot.savefig("beercan_fe55_voltage_scan.pdf")




####################################################################
#
# Americium data with the voltage scan
#
####################################################################
voltage_scan = {}
voltage_params = {}
files = glob.glob(os.path.abspath("../data/beer can/group 1/Beercan_Am*.mca"))
files.sort()
f, axes = pyplot.subplots(2,2)#, sharex=True, sharey=True)
i=0
j=0

for filepath in files: 
    fname = os.path.basename(filepath)
    voltage = float(fname.split("_")[2].replace("V",""))
    gain = float(fname.split("_")[-1].split(".")[0][1:])

    if gain < 20: continue

    # Get the data
    data, livetime = read_mca(filepath)
    #data = numpy.array(data[45:]) # Removing channels below 45 since that seems to efamctively remove noise

    #print numpy.sum(data/livetime)
    # Apply the calibration using the gain
    # Get this slope and intercept
    fitparams, fiterrors = calibration_params(gain, calib)
    
    channels = numpy.cumsum(numpy.ones_like(data))
    mvolts, charge_err = chantoV(channels, fitparams, fiterrors)

    voltage_scan[voltage] = [mvolts, charge_err, data/livetime, numpy.sqrt(data)/livetime]
    #fit, sigma = fitDoubleGaussians(mvolts, charge_err, data/livetime, numpy.sqrt(data)/livetime)

    #voltage_params[voltage] = [fit, sigma]
    axes[i,j].errorbar(x=mvolts, y=data/livetime, xerr=charge_err, yerr=numpy.sqrt(data)/livetime, label="Data")
    axes[i,j].set_xlim(numpy.min(mvolts), numpy.max(mvolts))
    axes[i,j].set_ylim(numpy.min(data/livetime), numpy.max(data/livetime))

    #fitvalues = double_gauss_fit(fit, mvolts)
    #axes[i,j].plot(mvolts, fitvalues, 'rx', label="Fit")
    
    i += 1
    if i == 2:
        j += 1
        i = 0

pyplot.tight_layout()
pyplot.xlabel("Charge (e)")
pyplot.ylabel("Events")
pyplot.savefig("beercan_am_voltage_scan.pdf")
