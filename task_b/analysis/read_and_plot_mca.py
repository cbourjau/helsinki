import numpy as np
import matplotlib.pyplot as plt
from utils import read_mca


if(__name__ == '__main__'):
    import sys
    file_names = sys.argv[1:]
    plt.figure()
    for f in file_names:
        data,live_time = read_mca(f)

        plt.step(range(len(data)),data,where='mid',label= "%s live time %f s"%(f,live_time))
    #channel np.array(range(len(data)))
    #mask = (channel>300) & (channel<400)
    #chargsum = sum(np.array(data)[mask])
    #print(chargsum)
    plt.legend()
    plt.xlabel("Channel",size=20)
    plt.ylabel("Count",size=20)
    plt.yscale('log')
    plt.show()
