#!/usr/bin/env python
# -*- coding: utf-8 -*-
'make distribution and K-S test of GRB redshift and luminosity '
__authou__='Wenjin Xie'

import numpy as np
import random
import  lumdist      #turn redshift(z) to Mpc
import randgenerator  
from math import *
import matplotlib.pyplot as plt
import pandas as pd
import time
import os
MPC=3.086*1e24    #Mpc--cm

#光度分布函数
def L_angel_fun(L0,k,j,v):      
    v1=v/180.0*np.pi
    if(j>v1):
        return L0
    else:
        L1=L0-k*(np.log10(v1)-np.log10(j))
        return L1
# 单位空间红移
def dVdz(z):     #dV/dz: c=3.0*1e10 cm  H0=7.1*1e6 cm/s/Mpc  omiga_m=0.7  omiga_ka=0.3   unit:Mpc3 to Gpc
    return (3.0*4*np.pi*lumdist.lumdist(z)**2)/(7.1*(1.0+z)**2*(0.3*(1.0+z)**3+0.7)**0.5)*1e-5
# 恒星形成率分布
def SFR2(zd): #  !Hasan 2008
    return 0.02*((1+zd)**(-3.4*10)+((1+zd)/5160.6)**(0.3*10)+((1+zd)/9.06)**(3.5*10))**(-1.0/10.0)
# 伽马暴爆发率与恒星形成率的关系
def R_z2(z):  #!Hasan 2008
    return SFR2(z)*(1.0+z)**1.5*dVdz(z)/(1.0+z)
# 角度分布函数
def angel_rand():
    return np.arccos(random.uniform(0.173,1.0))*180/np.pi   #0-80 degree
   # return np.arccos(random.uniform(0.766,1.0))*180/np.pi  #0-40 degree
   # return np.arccos(random.uniform(0.984,1.0))*180/np.pi  #0-10 degree
# k is structure index , core unit is degree , n is sample number, Lcore is core luminosity
def simulation(k=4.0,core=2.1,n=10000,Lcore=52.5):
    if(not(os.path.exists('final_sig_z2_figure'))):
        os.mkdir('final_sig_z2_figure')    
    if(not(os.path.exists('final_sig_z2_figure/png'))):
        os.mkdir('final_sig_z2_figure/png')     
  #  t0 = time.clock()
    Luminosity_rand=[]
    redshift_rand=[]
    vangle_rand = []  
    N = 10000
    n0= 0
    t0 = time.clock()
    while 1:
        L_rand = np.array([random.normalvariate(Lcore,0.5) for i in range(N)])
        v_angle = np.array([angel_rand() for i in range(N)])
        j_angle = core*3.1415/180.0
        L_obs =  [L_angel_fun(L_rand[i],k, j_angle, v_angle[i]) for i in range(len(v_angle))]
  #      L_obs =  [L_angel_fun(Lcore,k, j_angle, v_angle[i]) for i in range(len(v_angle))]
        z_rand= np.array(randgenerator.generator(R_z2,0.001,10.0,N))
        dl=[3.086*lumdist.lumdist(z_rand[i]) for i in range(len(z_rand))]   # 1Mpc=3.086E24  unit:1E24cm
        flux = [(L_obs[i]-48.0)-np.log10(3.0*4*np.pi*dl[i]**2) for i in range(len(L_obs))]  #k=3 roghtly
        for i in range(len(flux)):
            if flux[i] > random.uniform(-8,-7):  # instrument flux
                n0=n0+1 
                print (n0)
                Luminosity_rand.append(L_obs[i])
                redshift_rand.append(z_rand[i])
                vangle_rand.append(v_angle[i])
                if n0 >= n :
                    break                                                                 
        if n0 >= n:
            break
     #  write down the simulation data   
    t1 = time.clock()
    print ('calculated timd is', t1-t0)
    path='final_sig_z2_figure/'
    pathpng='final_sig_z2_figure/png/'
    datafile='all_data.xlsx'
    data=pd.read_excel(datafile)
    z=np.array(data['z'])
    Luminosity=np.array(data['Luminosity'])
    #---------------------------------------K-S test
    from scipy.stats import ks_2samp
    ks1,ks2=ks_2samp(Luminosity_rand,Luminosity)
    ks3,ks4=ks_2samp(redshift_rand,z)
  #  print (ks1,ks2)
 #  print (ks3,ks4)
    KS = ks2*ks4 
    
    with open(path+"k %5.2f c %5.2f L %5.2f KS %2.3f.txt"%(k,core,Lcore,KS),"w") as f:
        for i in range(len(redshift_rand)):
            i = str(redshift_rand[i])+'  '+str(Luminosity_rand[i])+'  '+str(vangle_rand[i])+'\n'
            f.writelines(i)   
    #------------------------plot fig
    redshift_rand=[np.log10(1.0+redshift_rand[i]) for i in range(len(redshift_rand))]
    z=[np.log10(1.0+z[i]) for i in range(len(z))]
    #print z
    # definitions for the axes
    left, width = 0.1, 0.45
    bottom, height = 0.1, 0.45
    bottom_h = left_h = left + width + 0.05
    
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.35]
    rect_histy = [left_h, bottom, 0.35, height]
    rect_histz = [bottom_h+0.05,bottom_h+0.05, 0.30, 0.30]
    # start with a rectangular Figure
    fig = plt.figure(1, figsize=(10, 10))
    axScatter = plt.axes(rect_scatter)
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)
    # the scatter plot:
    plt.axes(rect_scatter)
    plt.scatter(redshift_rand, Luminosity_rand,marker='o',c='',edgecolors='b',label='simulation data')
    plt.scatter(z, Luminosity,label='observer data',c='r')
    plt.xlabel('log(1+z)')
    plt.ylim(46.0,54.4)
    plt.xlim(0,1.0)
    plt.ylabel('Log L (erg/s)')
    plt.legend()
    
    plt.axes(rect_histy)
    weights3 = np.ones_like(Luminosity_rand)/float(len(Luminosity_rand))
    weights4 = np.ones_like(Luminosity)/float(len(Luminosity))
    plt.hist(Luminosity_rand, weights=weights3, bins=np.linspace(47,54,28), histtype='step',linestyle=('dashed'), normed=0, orientation='horizontal',label='simulation data')
    plt.hist(Luminosity, weights=weights4, bins=np.linspace(47,54,28), histtype='step', normed=0, orientation='horizontal',color='r',label='observer data')
    plt.ylabel('Log L (erg/s)')
    plt.xlabel(' probability')
    plt.title('Luminosity P$_{k-s}$=%3.2f'%ks2)
    plt.legend()
    
    plt.axes(rect_histx)
    weights1 = np.ones_like(redshift_rand)/float(len(redshift_rand))
    weights2 = np.ones_like(z)/float(len(z))
    plt.hist(redshift_rand, bins=np.linspace(0,1.0,20),weights=weights1,histtype='step', linestyle=('dashed'),normed=0, label='simulation data')
    plt.hist(z,bins=np.linspace(0,1.0,20),weights=weights2, color='r',histtype='step',normed=0, label='observer data')
    plt.xlabel(' log(1+z)')
    plt.ylabel(' probabilty')
    plt.title('Redshift P$_{k-s}$=%3.2f'%ks4)
    plt.legend()
    
    plt.axes(rect_histz)
    weights5 = np.ones_like(vangle_rand)/float(len(vangle_rand))
    plt.hist(vangle_rand, histtype='step',weights=weights5, normed=0, label='simulation angle')
    plt.ylabel(' probability')
    plt.xlabel('viewing angle(${\circ}$)')
    plt.legend()
    plt.savefig(pathpng+'k%5.3f c%5.3f L%2.3f KS%2.3f.png'%(k,core,Lcore,KS))
  #  plt.savefig(pathpng+'k%5.2f c%5.2f L%3.2f KS%2.2f.eps'%(k,core,Lcore,KS))
    plt.close()    
  #  t1 = time.clock()
 #   print ('calculated timd is', t1-t0)
    print ('k=%5.4f core=%5.4f  KS is %5.4f '%(k,core,KS))
    return k,core,KS,ks2,ks4,Lcore

#t0 = time.clock()
#simulation(k=6.0,core=0.50,n=100,Lcore=51.0)  
#t1 = time.clock()
#print ('calculated timd is', t1-t0)

if __name__=="__main__":
    simulation()