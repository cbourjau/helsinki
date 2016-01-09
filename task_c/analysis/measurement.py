import numpy as np
import utils
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#Gaussian function
def gauss_function(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

inpV = [10,20,30]
d,lt = utils.read_spe('../data/beata_spectrum1.Spe')
noise,noise_lt = utils.read_spe('../data/noise_spectrum1.Spe')
calib,trash = utils.read_spe('../data/calibration_right_2musec_10_20_30.Spe')
plt.step(range(len(calib)),calib,where='mid')
channels = np.arange(len(d))
#findings calibration peaks
peakest = 351
mask = (channels>=(peakest-5)) & (channels<=(peakest+5))
popt1, pcov = curve_fit(gauss_function, channels[mask], calib[mask],p0 = [10000,peakest, 1.0])
peakest = 700
mask = (channels>=(peakest-5)) & (channels<=(peakest+5))
popt2, pcov = curve_fit(gauss_function, channels[mask], calib[mask], p0 = [10000,peakest, 1.0])
peakest = 1051
mask = (channels>=(peakest-5)) & (channels<=(peakest+5))
popt3, pcov = curve_fit(gauss_function, channels[mask], calib[mask],p0 = [10000,peakest, 1.0])
peaks = [popt1[1],popt2[1],popt3[1]]
def calib_gen(inpV,peaks):
    def linear_fit(x,m,b): return m*x+b
    print(inpV)
    print(peaks)
    popt = np.polyfit(peaks,inpV,1)
    print(popt)
    def calib(channel):
        return popt[0]*channel+ popt[1]
    return calib
calibration = calib_gen(inpV,peaks)
v = calibration(channels)/10.07
plt.figure()
plt.plot(channels,calibration(channels))
plt.figure()
plt.step(v,noise/noise_lt,where='mid',label= "live time %f s"%(noise_lt))
plt.figure()
plt.step(v,d/lt-noise/noise_lt,where='mid',label= "live time %f s"%(lt))
#plt.legend()
plt.xlabel("Charge equivalent of a MIP",size=20)
plt.ylabel("Rate [Hz]",size=20)
plt.yscale('log')
plt.xlim((0,1.5))
plt.show()
