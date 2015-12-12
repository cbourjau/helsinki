import numpy as np
import matplotlib.pyplot as plt
from utils import read_mca

data20,live_time = read_mca("../data/calib_detector_gain20_5peaks.mca")
data50,live_time = read_mca("../data/calib_detector_gain50_5peaks.mca")
data100,live_time = read_mca("../data/calib_detector_gain100_4peaks.mca")
channels = np.arange(len(data20))
#Pulse heights from Christians notes
pulseHeights20 = np.array([10.3,20.3,39.6,70.,110.])#mV
pulseHeights50 = np.array([5.,10.0,20.,40.,50.])#mV
pulseHeights100 = np.array([2.,5.0,10.0,19.6])#mV
#By eye estimates where the peaks roughly are
peaksrough20 = np.array([25,72,152,276,436])
peaksrough50 = np.array([36,90,193,388,494])
peaksrough100 = np.array([15,76,175,371])

plt.figure(1)
plt.step(channels,data20,where='mid',label= " Gain 20")
plt.figure(2)
plt.step(channels,data50,where='mid',label= " Gain 50")
plt.figure(3)
plt.step(channels,data100,where='mid',label= "Gain 100 ")

plt.legend()
plt.xlabel("Channel",size=20)
plt.ylabel("Count",size=20)
plt.yscale('log')

#Gaussian function
def gauss_function(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))


# program
from scipy.optimize import curve_fit
#Finding peaks for gain 20
plt.figure(1)
peaks20 = list()
plt.title("Gain 20")
for i in range(len(peaksrough20)):
    peakest = peaksrough20[i]
    mask = (channels>=(peakest-10)) & (channels<=(peakest+10))
    popt, pcov = curve_fit(gauss_function, channels[mask], data20[mask],p0 = [1000,peakest, 1.0])
    #print(popt)
    peaks20.append(popt)
    x = np.linspace(np.min(channels[mask]),np.max(channels[mask]),100)
    plt.plot(x,gauss_function(x,popt[0],popt[1],popt[2]),color='red',linestyle = '--',lw=2)
plt.ylim((1,1e5))

#Finding peaks for gain 50
plt.figure(2)
peaks50 = list()
plt.title("Gain 50")
for i in range(len(peaksrough50)):
    peakest = peaksrough50[i]
    mask = (channels>=(peakest-15)) & (channels<=(peakest+15))
    popt, pcov = curve_fit(gauss_function, channels[mask], data50[mask],p0 = [1000,peakest, 2.0])
    #print(popt)
    peaks50.append(popt)
    x = np.linspace(np.min(channels[mask]),np.max(channels[mask]),100)
    plt.plot(x,gauss_function(x,popt[0],popt[1],popt[2]),color='red',linestyle = '--',lw=2)

plt.ylim((1,1e5))

#Finding peaks for gain 100
plt.figure(3)
peaks100 = list()
plt.title("Gain 100")
for i in range(len(peaksrough100)):
    peakest = peaksrough100[i]
    if(i==0):
        low = 5
    else:
        low = 15
    mask = (channels>=(peakest-low)) & (channels<=(peakest+20))
    popt, pcov = curve_fit(gauss_function, channels[mask], data100[mask],p0 = [1000,peakest, 2.0])
    #print(popt)
    peaks100.append(popt)
    x = np.linspace(np.min(channels[mask]),np.max(channels[mask]),100)
    plt.plot(x,gauss_function(x,popt[0],popt[1],popt[2]),color='red',linestyle = '--',lw=2)

plt.ylim((1,1e5))


peaks20 = zip(*peaks20)
peaks50 = zip(*peaks50)
peaks100 = zip(*peaks100)

def lin_fun(x,k,m):
    return x*k+m
#Fitting a line to the input voltage and channel peak
popt20, pcov20 = curve_fit(lin_fun, pulseHeights20, peaks20[1] ,p0 = [10.0,0])
popt50, pcov50 = curve_fit(lin_fun, pulseHeights50, peaks50[1] ,p0 = [10.0,0])
popt100, pcov100 = curve_fit(lin_fun, pulseHeights100, peaks100[1] ,p0 = [10.0,0])
plt.figure()
x = np.linspace(0,max(pulseHeights20),100)
plt.errorbar(pulseHeights20,peaks20[1],yerr=peaks20[2],color = 'blue',label = "Gain 20, channel/mV %f"%popt20[0],ls = 'None')
plt.plot(x,lin_fun(x,popt20[0],popt20[1]),color = 'blue')

x = np.linspace(0,max(pulseHeights50),100)
plt.errorbar(pulseHeights50,peaks50[1],yerr=peaks50[2],color = 'red',label = "Gain 50, channel/mV %f"%popt50[0],ls = 'None')
plt.plot(x,lin_fun(x,popt50[0],popt50[1]),color = 'red')

x = np.linspace(0,max(pulseHeights100),100)
plt.errorbar(pulseHeights100,peaks100[1],yerr=peaks100[2],color = 'green',label = "Gain 100, channel/mV %f"%popt100[0],ls = 'None')
plt.plot(x,lin_fun(x,popt100[0],popt100[1]),color = 'green')


plt.ylabel("Channel")
plt.xlabel("Input pulse height [mV]")
plt.legend()

#Fitting a line to the channel peak and input voltage (inverse of what was done above)
popt20, pcov20 = curve_fit(lin_fun, np.array(peaks20[1]), np.array(pulseHeights20), p0 = [1/10.,0])
popt50, pcov50 = curve_fit(lin_fun, np.array(peaks50[1]), np.array(pulseHeights50) ,p0 = [1/10.0,0])
popt100, pcov100 = curve_fit(lin_fun, np.array(peaks100[1]), np.array(pulseHeights100) ,p0 = [1/10.0,0])
plt.figure()
x = np.linspace(0,512,100)
plt.errorbar(peaks20[1],pulseHeights20,xerr=peaks20[2],color = 'blue',label = "Gain 20, mV/channel %f"%popt20[0],ls = 'None')
plt.plot(x,lin_fun(x,popt20[0],popt20[1]),color = 'blue')

plt.errorbar(peaks50[1],pulseHeights50,xerr=peaks50[2],color = 'red',label = "Gain 50, mV/channel %f"%popt50[0],ls = 'None')
plt.plot(x,lin_fun(x,popt50[0],popt50[1]),color = 'red')

plt.errorbar(peaks100[1],pulseHeights100,xerr=peaks100[2],color = 'green',label = "Gain 100, mV/channel %f"%popt100[0],ls = 'None')
plt.plot(x,lin_fun(x,popt100[0],popt100[1]),color = 'green')


plt.xlabel("Channel")
plt.ylabel("Input pulse height [mV]")
plt.legend()
#Saving the calibraton in a pickle file
import pickle
calibration = dict()
calibration['gain20'] = popt20
calibration['gain50'] = popt50
calibration['gain100'] = popt100

f = open("../data/calibration.pkl",'w')
pickle.dump(calibration,f)

plt.show()
