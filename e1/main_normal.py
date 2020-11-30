# -*- coding: utf-8 -*-
"""
多线程运行脚本
"""
import final_sigle_z2_figure_normalL as e1
import  numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri


def proc(N=30):
  #  L1 = []
    print('Run task %s (%s)...' % (N, os.getpid()))
    for i in range(N):
        Lc = np.random.uniform(51,54)
        k = np.random.uniform(2,6.0)
        core = np.random.uniform(0.5,20)
        print ('k=%5.2f,core=%5.2f,count=%5.0f'%(k,core,i))
        L = e1.simulation(k,core,n=1500,Lcore=Lc)  #retrun k,core,KS,ks_L,ks_z
        with open("randL_ks_data.txt","a") as f:
            i = str(L[0])+'  '+str(L[1])+'  '+str(L[5])+'  '+str(L[2])+'  '+str(L[3])\
            +'  '+str(L[4])+'\n'
            f.writelines(i)
  #  print('Task %s runs %0.2f seconds.' % (N, (end - start)))
     #   L1.append(L)
from multiprocessing import Pool
import multiprocessing
import os
import time
if __name__ == '__main__':
    starttime = time.time()
    processes = []
    for i in range(0,30):
        p = multiprocessing.Process(target=proc, args=(666,))
        processes.append(p)
        p.start()
    for process in processes:
        process.join()
    print('That took {} seconds'.format(time.time() - starttime))
    
    