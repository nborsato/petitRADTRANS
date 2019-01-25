from __future__ import print_function

import numpy as np
import struct

def read_bin(filename):

    # open the file
    file = open(filename,mode='rb')
    # read the content
    cont = file.read()

    '''
    Every row is delimited by an integer number.
    If we have two columns containing two doubles,
    then the number of bytes per row is 4+8+8+4 = 24
    '''
    # get the number of lines
    points = len(cont)/24
    #points = len(cont)/16
    # create arrays of the appropriate length
    x = np.ones(int(points)) # first column
    y = np.ones(int(points)) # second column
    
    # read the binary data into the arrays
    for i in range(int(points)):
        #test = struct.unpack('iddi',cont[i*24:(i+1)*24])
        test = struct.unpack('dd',cont[i*24+4:(i+1)*24-4])
        x[i],y[i] = test[0],test[1]
        #print(test)

    file.close()
        
    return x,y
