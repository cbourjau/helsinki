import numpy as np
import matplotlib.pyplot as plt
import pickle
from scipy.optimize import curve_fit
from utils import *

calibration = pickle.load(open("../data/calibration.pkl",'r'))

g = np.linspace(0,200,100)

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
