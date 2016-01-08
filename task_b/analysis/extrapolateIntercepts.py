import numpy as np
import matplotlib.pyplot as plt
import pickle
from scipy.optimize import curve_fit
from utils import *
from mytools import plot_setup
import matplotlib
plot_setup.setup_mpl(matplotlib)
calibration = pickle.load(open("../data/calibration.pkl",'r'))

g = np.linspace(0,200,100)

chans = np.linspace(0,512,10000)
plt.figure(1,figsize = (12,8))
plt.figure(2,figsize = (12,8))
for g in [20,50,100,200,500]:
    p ,p_err = calibration_params(g,calibration)
    print(g)
    print("&$%f \pm %g$ & $%g \pm %g$"%(p[0],p_err[0],p[1],p_err[1]))
    V,V_err = chantoE(chans,p,p_err)
    plt.figure(1)
    plt.errorbar(chans,V,yerr=V_err,label="Gain %d"%g)
    plt.figure(2)
    plt.plot(chans,V_err/V,label="Gain %d"%g)
plt.figure(1)
plt.ylabel("Number of electrons",size=20)
plt.xlabel("Channel",size=20)
plt.legend(loc='best')
plt.savefig("CalibrationErrorDense.png")
plt.figure(2)
plt.ylabel("Relative error",size=20)
plt.xlabel("Channel",size=20)
plt.legend(loc='best')
plt.yscale('log')
plt.savefig("CalibrationErrorRelative.png")

plt.show()
