import numpy as np

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
            data.append(float(lines[i]))
    return  np.array(data), live_time
