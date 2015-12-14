import numpy as np
import matplotlib.pyplot as plt
import pickle
from scipy.optimize import curve_fit
from utils import *

calibration = pickle.load(open("../data/calibration.pkl",'r'))

slopes = list()
intercepts = list()
slopes_err = list()
intercepts_err = list()

for k,v in calibration.iteritems():
    slopes.append(v[0][0])
    intercepts.append(v[0][1])
    slopes_err.append(np.sqrt(v[1][0][0]))
    intercepts_err.append(np.sqrt(v[1][1][1]))
gains = [20,50,100]

plt.figure()
plt.errorbar(gains,slopes,yerr=slopes_err)
plt.errorbar(gains,intercepts,yerr=intercepts_err)
plt.xlim((0,150))



g = np.linspace(min(gains),200,100)

pop, pcov = curve_fit(pow, gains, slopes, p0 = [1.,-1.0], sigma = slopes_err)
plt.plot(g,pow(g,pop[0],pop[1]))

chans = np.linspace(0,512,100)
plt.figure()
for g in [20,50,100,200,500]:
    p ,p_err = calibration_params(g,calibration)
    V,V_err = chantoE(chans,p,p_err)
    plt.errorbar(chans,V,yerr=V_err,label="Gain %d"%g)
plt.ylabel("Number of electrons")
plt.xlabel("Channel")
plt.legend()
plt.show()
