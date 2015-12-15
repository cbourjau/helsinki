import numpy as np
import matplotlib.pyplot as plt
from utils import read_mca
figsize = (12,8)
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
plt.step(channels,data20, lw = 2, where='mid',label= "Calibration Data")
plt.figure(2)
plt.step(channels,data50, lw = 2, where='mid',label= "Calibration Data")
plt.figure(3)
plt.step(channels,data100, lw = 2, where='mid',label= "Calibration Data")


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
plt.plot([],[],color='red',linestyle = '--',lw=2,label='Gaussain fit to peaks')

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
plt.plot([0],[0],color='red',linestyle = '--',lw=2,label='Gaussain fit to peaks')

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
plt.plot([],[],color='red',linestyle = '--',lw=2,label='Gaussain fit to peaks')


for i in range(1,4):
    plt.figure(i)
    plt.xlabel("Channel",size=20)
    plt.ylabel("Count",size=20)
    plt.yscale('log')
    plt.legend(loc='best')
    plt.ylim((100,1e5))
    plt.savefig("CalibratioPeaks%d.png"%i)

peaks20 = zip(*peaks20)
peaks50 = zip(*peaks50)
peaks100 = zip(*peaks100)

def lin_fun(x,k,m):
    return x*k+m


#Fitting a line to the input voltage and channel peak
popt20, pcov20 = curve_fit(lin_fun, pulseHeights20, peaks20[1] ,p0 = [10.0,0])


popt50, pcov50 = curve_fit(lin_fun, pulseHeights50, peaks50[1] ,p0 = [10.0,0])
popt100, pcov100 = curve_fit(lin_fun, pulseHeights100, peaks100[1] ,p0 = [10.0,0])
plt.figure(figsize=figsize)
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
plt.legend(loc='best')

#Fitting a line to the channel peak and input voltage (inverse of what was done above)
pulseheight_er = 0.3
popt20, pcov20 = curve_fit(lin_fun, np.array(peaks20[1]), np.array(pulseHeights20), p0 = [1/10.,0],sigma = np.ones(5)*pulseheight_er)
popt50, pcov50 = curve_fit(lin_fun, np.array(peaks50[1]), np.array(pulseHeights50), p0 = [1/10.0,0],sigma = np.ones(5)*pulseheight_er)
popt100, pcov100 = curve_fit(lin_fun, np.array(peaks100[1]), np.array(pulseHeights100), p0 = [1/10.0,0],sigma = np.ones(4)*pulseheight_er)

plt.figure(figsize=figsize)
x = np.linspace(0,512,100)
plt.errorbar(peaks20[1],pulseHeights20,yerr=np.ones(5)*pulseheight_er,color = 'blue',label = "Gain 20, mV/channel %f"%popt20[0],ls = 'None')
plt.plot(x,lin_fun(x,popt20[0],popt20[1]),color = 'blue')

plt.errorbar(peaks50[1],pulseHeights50,yerr=np.ones(5)*pulseheight_er,color = 'red',label = "Gain 50, mV/channel %f"%popt50[0],ls = 'None')
plt.plot(x,lin_fun(x,popt50[0],popt50[1]),color = 'red')

plt.errorbar(peaks100[1],pulseHeights100,yerr=np.ones(4)*pulseheight_er,color = 'green',label = "Gain 100, mV/channel %f"%popt100[0],ls = 'None')
plt.plot(x,lin_fun(x,popt100[0],popt100[1]),color = 'green')


plt.xlabel("Channel")
plt.ylabel("Input pulse height [mV]")
plt.legend(loc='best')
plt.savefig("CalibrationFit.png")

#Fitting the slopes and Intercepts to the gains for extrapolation and interpolation
gains = np.array([20,50,100])
k = np.array([popt20[0],popt50[0],popt100[0]])
m = np.array([popt20[1],popt50[1],popt100[1]])

#Errors taken from the fits
k_err =  np.sqrt(np.array([pcov20[0][0],pcov50[0][0],pcov100[0][0]]))
m_err =  np.sqrt(np.array([pcov20[1][1],pcov50[1][1],pcov100[1][1]]))

def powe(x,n):
    return n*x**(-1.0)
def poweg(x,n,i):
    return n*x**(i)
def powegdi(x,n,i):
    return n*np.log(x)*x**i
def powegdn(x,n,i):
        return x**i
poptgain3, pcovgain3 =curve_fit(poweg, gains, k , sigma= k_err,p0 = [1.0,-1])


plt.figure(figsize=figsize)
plt.errorbar(gains,k,yerr=k_err,ls='None',marker='o')
g = np.linspace(1,500,1000)
plt.plot(g,poweg(g,poptgain3[0],poptgain3[1]),color='blue',
        label='normalization = %f'%poptgain3[0]+r"$\pm$ "+"%f, index = %f"%(np.sqrt(pcovgain3[0][0]),poptgain3[1]) +r"$\pm$ "+"%f"%np.sqrt(pcovgain3[1][1]))

plt.yscale('log')
plt.xlabel("gain")
plt.ylabel("mV/channel")
plt.legend(loc='best')



plt.figure(figsize=figsize)
plt.title("Intercepts")
poptgainint, pcovgainint =curve_fit(poweg, gains, m , sigma = m_err, p0 = [1.0,-1])
plt.plot(g,poweg(g,poptgainint[0],poptgainint[1]),color='blue',
        label='normalization = %f'%poptgainint[0]+r"$\pm$ "+"%f, index = %f"%(np.sqrt(pcovgainint[0][0]),poptgainint[1]) +r"$\pm$ "+"%f"%np.sqrt(pcovgainint[1][1]))
plt.errorbar(gains,m,yerr=k_err,ls='None',marker='o')
plt.xlabel("gain")
plt.ylabel("intercept [mV]")
plt.legend(loc='best')


def calibration_params(gain, slopes, slope_er, intercept, intercept_er):
    err_slope = np.sqrt(powegdi(gain,slopes[0],slopes[1])**2*slope_er[0][0] +
                        powegdn(gain,slopes[0],slopes[1])**2*slope_er[1][1] +
                        2*powegdn(gain,slopes[0],slopes[1])*powegdi(gain,slopes[0],slopes[1])*slope_er[1][0])
    err_intercept = np.sqrt( powegdi(gain,intercept[0], intercept[1])**2*intercept_er[0][0]+
                             powegdn(gain,intercept[0], intercept[1])**2*intercept_er[1][1]+
                             2*powegdn(gain,intercept[0],intercept[1])*powegdi(gain,intercept[0],intercept[1])*intercept_er[1][0])
    return (poweg(gain,slopes[0],slopes[1]),poweg(gain,intercept[0],intercept[1]) ),(err_slope,err_intercept)

#Saving the calibraton in a pickle file
import pickle
calibration = dict()

calibration["gain20"] = (popt20,pcov20)
calibration["gain50"] = (popt50,pcov50)
calibration["gain100"] = (popt100,pcov100)

f = open("../data/calibration.pkl",'w')
pickle.dump(calibration,f)

print("Parameters for gain 5 ",calibration_params(5,poptgain3, pcovgain3,poptgainint, pcovgainint))
print("Parameters for gain 200 ",calibration_params(200,poptgain3, pcovgain3,poptgainint, pcovgainint))
print("Parameters for gain 500 ",calibration_params(500,poptgain3, pcovgain3,poptgainint, pcovgainint))

def chantoV(chan, param, er_param):
    err = np.sqrt(chan**2*er_param[0]**2+er_param[1]**2)
    return (lin_fun(chan,param[0],param[1]),err)

chan = np.linspace(0,512,100)
plt.figure(figsize=figsize)
for g in [5,20,100,200,500]:
    p,p_err = calibration_params(g,poptgain3, pcovgain3,poptgainint, pcovgainint)
    f,err =chantoV(chan,p,p_err)
    plt.errorbar(chan,f,err,label="Gain %d"%g)
plt.legend(loc='best')
plt.xlabel("channel")
plt.ylabel("infered Voltage [mV]")

plt.figure(figsize=figsize)
perr =  [np.sqrt(pcov20)[0][0],np.sqrt(pcov20)[1][1]]
f,err =chantoV(chan,calibration['gain20'][0],perr)
plt.errorbar(chan,f,err,label="Gain 20")
perr =  [np.sqrt(pcov50)[0][0],np.sqrt(pcov50)[1][1]]
f,err =chantoV(chan,calibration['gain50'][0],perr)
plt.errorbar(chan,f,err,label="Gain 50")


perr =  [np.sqrt(pcov100)[0][0],np.sqrt(pcov100)[1][1]]
f,err =chantoV(chan,calibration['gain100'][0],perr)
plt.errorbar(chan,f,err,label="Gain 100")
plt.xlabel("channel")
plt.ylabel("infered Voltage [mV]")
plt.legend(loc='best')
plt.savefig("CalibrationFunctionError.png")

plt.show()
