import numpy as np


####################################################################
# Read the data from the .spe file and return less useless stuff
####################################################################
def read_spe(filename):
    '''A simple function to read .spe files
        - filename    path/filename to the .mca files
        returns
        data and live time
    '''
    with open(filename,'r') as f:
        lines = f.readlines()

        index = 0

        live_time = float(lines[9].split()[0])
        bindex = 12+int(lines[11].split()[1])
        data = list()
        for i in range(12,bindex):
            data.append(float(lines[i]))
    return  np.array(data), live_time
