#!/usr/bin/env python
# -*- coding: utf-8 -*-

'a generator than use Monte Calo method to \
    generete romdom by density distribution'

__authou__='Wenjin Xie'

import numpy as np
import numpy.random as random
import matplotlib.pyplot as plt

# 找出方程里面的最大值
def findmax(func,a,b):   
    x=np.arange(a,b,(b-a)/10000.0)
    y=[func(x[i]) for i in range(len(x))]
    return max(y) 

# 舍取法产生对应方程随机数
def generator(func,a,b,N=1000,*args):    
    A=1.01*findmax(func, a, b)
    def rand():
        while True:
            X = random.uniform(a,b)
            Y = random.uniform(0.0,1.0)*A
            Z = func(X)
            if Y<Z:
                return X 
    rand_data=[rand() for i in range(N)]
    return rand_data
   # return np.array(rand_data).reshape(-1,1)

def plotrand(a,bin_num,name=None,yscale=None,xscale=None):
    import matplotlib.pyplot
    plt.hist(a,bins=bin_num,normed=1)
    if yscale is None:
        pass
    else:
        plt.yscale('%s'%yscale)
    if xscale is None:
        pass
    else:
        plt.xscale('%s'%xscale)    
    if name==None:
        pass
    else:
        plt.title('%s'%name)
        plt.savefig("%s.png"%name)
    plt.show()
    plt.close()

if __name__=="__main__":
    generator(func, a, b)