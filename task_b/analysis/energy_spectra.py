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

lines = [5.19e3, 5.9e3] # keV

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
    
    print gain
    if gain < 20: continue

    #fitparams, fiterrors = calibration_params(gain, calib)
    #charges, charge_err = chantoE(data, fitparams, fiterrors)

    # Save the noise as channel number, error, rate, error
    noise_scan[gain] = [numpy.cumsum(numpy.ones_like(data)), data/livetime, numpy.sqrt(data)/livetime]

    
#noise_fits, cuts = fitNoise(noise_scan)

# Plot the spectra
hist1, bins, junk = pyplot.hist(noise_scan[50][0], weights=noise_scan[50][2]/numpy.sum(noise_scan[50][2]),
                              bins = numpy.cumsum(numpy.ones(500)),
                              label="Gain = %i" % 50, color='r', histtype='step')
hist2, bins, junk = pyplot.hist(noise_scan[200][0], weights=noise_scan[200][2]/numpy.sum(noise_scan[200][2]),
                               bins = numpy.cumsum(numpy.ones(500)),
                               label="Gain = %i" % 200, color='g', histtype='step')
hist3, bins, junk = pyplot.hist(noise_scan[500][0], weights=noise_scan[500][2]/numpy.sum(noise_scan[500][2]), 
                               bins = numpy.cumsum(numpy.ones(500)),
                               label="Gain = %i" % 500, color='b', histtype='step')

pyplot.figure()
pyplot.hist(bins[:-1], bins=bins, weights=[numpy.sum(hist1[:i]) for i in range(len(hist1))], color='r', alpha=0.6, linewidth=3, histtype='step', label="Gain = %i" % 50)

pyplot.hist(bins[:-1], bins=bins, weights=[numpy.sum(hist2[:i]) for i in range(len(hist2))], color='g', alpha=0.6, linewidth=3, histtype='step', label="Gain = %i" % 100)

pyplot.hist(bins[:-1], bins=bins, weights=[numpy.sum(hist3[:i]) for i in range(len(hist3))], color='b', alpha=0.6, linewidth=3, histtype='step', label="Gain = %i" % 500)

pyplot.legend(loc='lower right')
pyplot.grid()
pyplot.xlim(0,500)
#pyplot.yscale('log')
pyplot.xlim(min(noise_scan[gain][0]), max(noise_scan[gain][0]))
#pyplot.ylim(min(noise_scan[gain][2]), max(noise_scan[gain][2]))
#pyplot.xlim(0, 75)
pyplot.xlabel("Channel Number")
pyplot.ylabel("Cumulative Fraction of Events")
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
f, axes = pyplot.subplots(3,3)#, sharex=True, sharey=True)
i=0
j=0

for filepath in files: 
    fname = os.path.basename(filepath)
    voltage = float(fname.split("_")[2].replace("HV",""))
    gain = float(fname.split("_")[-1][4:-4])

    if gain < 20: continue

    # Get the data
    data, livetime = read_mca(filepath)
    data = numpy.array(data[45:]) # Removing channels below 45 since that seems to effectively remove noise
    print livetime

    # Apply the calibration using the gain
    # Get this slope and intercept
    fitparams, fiterrors = calibration_params(gain, calib)
    
    channels = numpy.cumsum(numpy.ones_like(data))
    mvolts, charge_err = chantoV(channels, fitparams, fiterrors)

    voltage_scan[voltage] = [mvolts, charge_err, data/livetime, numpy.sqrt(data)/livetime]
    fit, sigma = fitDoubleGaussians(mvolts, charge_err, data/livetime, numpy.sqrt(data)/livetime)

    data = numpy.array(data[45:]) # Removing channels below 45 since that seems to effectively remove noise

    voltage_params[voltage] = [fit, sigma]
    #axes[i,j].errorbar(x=mvolts, y=data/livetime, xerr=charge_err, yerr=numpy.sqrt(data)/livetime, label="Data")
    #axes[i,j].set_xlim(numpy.min(mvolts), numpy.max(mvolts))
    #axes[i,j].set_ylim(numpy.min(data/livetime), numpy.max(data/livetime))

    voltage_scan[voltage].append(double_gauss_fit(fit, mvolts))
    axes[i,j].plot(mvolts, double_gauss_fit(fit, mvolts), 'rx', label="Fit")
    
    i += 1
    if i == 3:
        j += 1
        i = 0

pyplot.tight_layout()
pyplot.xlabel("Charge (e)")
pyplot.ylabel("Events")
pyplot.savefig("voltage_scan.pdf")


pyplot.figure()
bins = numpy.linspace(-2,2, 40, False)
for v in voltage_scan.keys():
    cut = voltage_scan[v][2] != 0
    expectation = voltage_scan[v][-1]
    print numpy.sum(voltage_scan[v][-1][cut])
    print numpy.sum(voltage_scan[v][2][cut])
    pyplot.hist((voltage_scan[v][2][cut] - expectation[cut])/expectation[cut],
                bins=bins, alpha=0.6, histtype='stepfilled', label=str(v)+" kV")
pyplot.xlabel("(Measured Events - Expected Events)/Expected Events")
pyplot.ylabel("Rate (Hz)")
pyplot.legend(loc='upper left')
pyplot.grid()
pyplot.savefig("DoubleGaussianErrors.pdf")



# Now we have values. Hurrah!
# lets plot the peak positions vs the voltage just for fun
voltages = voltage_params.keys()
voltages.sort()
peak1 = numpy.array([voltage_params[v][0][1] for v in voltages])
peak2 = numpy.array([voltage_params[v][0][4] for v in voltages])
sigma1 = numpy.array([voltage_params[v][0][2] for v in voltages])
sigma2 = numpy.array([voltage_params[v][0][5] for v in voltages])

pyplot.figure()
pyplot.plot(voltages, peak2, color='r', label="Fe55, 5.19 keV")
pyplot.fill_between(voltages, peak1-sigma1, peak1+sigma1, color='r', alpha=0.4)
pyplot.plot(voltages, peak1, color='b', label="Fe55, 5.90 keV")
pyplot.fill_between(voltages, peak2-sigma2, peak2+sigma2, color='b', alpha=0.4)
pyplot.legend(loc='upper left')
pyplot.yscale("log")
pyplot.grid()
pyplot.xlabel("High Voltage (kV)")
pyplot.ylabel("Peak Position (eV)")
pyplot.savefig("peak_positions.pdf")

# And do the resolution (ie, sigma) vs voltage
# Now convert! always assume the higher energy line is "correct"
[5.19e3, 5.9e3]

delta = peak1 - peak2
deltaE = 5.9e3 - 5.19e3
const_shift = 5.9e3 - numpy.multiply(deltaE/delta,peak1)

highvoltages = voltage_scan.keys()
#sortindices  = numpy.argsort(highvoltages)
#highvoltages = highvoltages[sortindices]

e1 = numpy.multiply(peak1, deltaE/delta) + const_shift
s1 = numpy.multiply(sigma1, deltaE/delta)
e2 = numpy.multiply(peak2, deltaE/delta) + const_shift
s2 = numpy.multiply(sigma2, deltaE/delta)

#e1 = e1[sortindices]
#s1 = s1[sortindices]
#e2 = e2[sortindices]
#s2 = s2[sortindices]

pyplot.figure()
f, axes = pyplot.subplots(1,2)
pyplot.figure()
pyplot.errorbar(e2/1e3, highvoltages, xerr=s2/1e3, color='r', label="Fe55, 5.19 keV", linewidth=3)
pyplot.errorbar(e1/1e3, highvoltages, xerr=s1/1e3, color='b', label="Fe55, 5.90 keV", linewidth=3)
#pyplot.fill_betweenx(highvoltages, e2-s2, e2+s2, color='r')
#pyplot.fill_betweenx(highvoltages, e1-s1, e1+s1, color='b')
pyplot.legend(loc='upper center')
#pyplot.yscale("log")
pyplot.grid()
pyplot.ylabel("High Voltage (kV)")
pyplot.xlabel("Peak Position (keV)")
pyplot.savefig("energy_positions.pdf")


####################################################################
#
# Cesium analysis
# These values are mV/channel and are for each gain in the form
# of a {gain: [slope, intercept]} from a plot of the mV vs channel
#
####################################################################
data, livetime = read_mca(os.path.abspath("../data/cesiumHV0.96_gain50.mca"))
print livetime
data = data[45:]
print numpy.sum(data)
gain = 50
high_voltage = 0.96

delta = peak1[highvoltages == high_voltage] - peak2[highvoltages == high_voltage]
deltaE = 5.9e3 - 5.19e3
const_shift = 5.9e3 - numpy.multiply(deltaE/delta,peak1)

channels = numpy.cumsum(numpy.ones_like(data))
fitparams, fiterrors = calibration_params(gain, calib)
mvolts, charge_err = chantoV(channels, fitparams, fiterrors)

print delta, deltaE, const_shift
energy = numpy.multiply(mvolts, deltaE/delta) + const_shift[highvoltages == high_voltage]

pyplot.figure()
pyplot.scatter(energy/1e3, data/(livetime))
pyplot.grid()
pyplot.xlim(min(energy/1e3), max(energy/1e3))
pyplot.xlabel("Energy (keV)")
pyplot.ylabel("Rate (Hz)")
pyplot.savefig("cesium.pdf")
