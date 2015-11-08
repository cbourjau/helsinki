import numpy as np
import matplotlib.pyplot as plt


def read_mca(filename):
    '''A simple function to read .mca files
        - filename    path/filename to the .mca files
        returns
        data and live time
    '''
    with open(filename,'r') as f:
        lines = f.readlines()

        index = 0
        for l in lines:
            if("LIVE_TIME" in l):
                live_time = float(l.split()[2])
            if("<<DATA>>" in l):
                break
            index += 1
        data = list()
        for i in range(index+1,len(lines)-1):
            data.append(int(lines[i]))
    return data, live_time

if(__name__ == '__main__'):
    import sys
    file_names = sys.argv[1:]
    plt.figure()
    for f in file_names:
        data,live_time = read_mca(f)

        plt.step(range(len(data)),data,where='mid',label= f)
    plt.legend()
    plt.xlabel("Channel",size=20)
    plt.ylabel("Count",size=20)
    plt.show()
