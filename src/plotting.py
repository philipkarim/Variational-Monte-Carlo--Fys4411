# Common imports
import os
import numpy as np
import pandas as pd
from pandas import DataFrame
import matplotlib.pyplot as plt
#from scipy import *


def data_path(DATA_ID, dat_id):
    return os.path.join(DATA_ID, dat_id)

def variationalpha():
    #Filenames
    fn=['07N10Dim3.txt', '03N10Dim3.txt', '07N100Dim3.txt', '03N100Dim3.txt']
    fn_interact=['N10Dim3']
    #Folders
    folder_noninteract = ["Results/GDalpha/noninteract/bruteforce/analytic"]
    folder_interact = ["Results/GDalpha/interact/bruteforce/numeric"]

    #infile = open(data_path(folder_noninteract[0], fn[4]),'r')

    infile = np.loadtxt(data_path(folder_noninteract[0], fn[1]))

    #print(infile[:,1])


    #plt.plot(infile[:,0], infile[:,1], 'ro')
    #plt.show()

    """
    Plotting variation of alpha
    """
    #Dim=3, particles=10, noninteracting, 2^18 steps
    a_nonint=  [0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7]
    el_nonint= [16.98, 15.92, 15.38, 15.09, 15, 15.04, 15.26, 15.47, 15.92]

    plt.plot(a_nonint, el_nonint)
    plt.plot(a_nonint, el_nonint, 'ro', markersize=3)
    plt.xlabel(r'$\alpha$',fontsize=16)
    plt.ylabel(r'$E_L (\hbar \omega)$',fontsize=16)
    plt.show()

    #Dim=3, particles=10, interacting, 2^16 steps
    a_int=  [0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7]
    el_int= [28.13, 25.76, 24.83, 24.55, 24.39, 24.54, 24.72, 25.15, 25.79]

    plt.plot(a_int, el_int)
    plt.plot(a_int, el_int, 'ro', markersize=3)
    plt.xlabel(r'$\alpha$',fontsize=16)
    plt.ylabel(r'$E_L (\hbar \omega)$',fontsize=16)
    plt.show()

    return


def plottsteps():
    #Filenames
    fn_bf=['steplength0.1N10Dim3.txt', 'steplength0.500000N10Dim3.txt', 'steplength1N10Dim3.txt']
    fn_is=['timestep0.500000N10Dim3.txt', 'timestep0.500000N10Dim3.txt', 'timestep0.500000N10Dim3.txt', 'timestep0.500000N10Dim3.txt']

    #Folders
    folder= ["Results/steps/bruteforce/", "Results/GDalpha/interact/bruteforce/numeric"]

    #infile = open(data_path(folder_noninteract[0], fn[4]),'r')

    infile = np.loadtxt(data_path(folder[0], fn_bf[0]))
    infile2 = np.loadtxt(data_path(folder[0], fn_bf[1]))
    infile3 = np.loadtxt(data_path(folder[0], fn_bf[2]))

    #print(infile[:,1])


    plt.plot(infile[:,0], infile[:,1], label='Step Length=0.1')
    plt.plot(infile2[:,0], infile2[:,1], label='Step Length=0.5')
    plt.plot(infile3[:,0], infile3[:,1], label='Step Length=1')
    plt.xlabel('Monte Carlo steps',fontsize=16)
    plt.ylabel(r'$E_L (\hbar \omega)$',fontsize=14)
    plt.show()


    infile4 = np.loadtxt(data_path(folder[1], fn_bf[0]))
    infile5 = np.loadtxt(data_path(folder[1], fn_bf[1]))
    infile6 = np.loadtxt(data_path(folder[1], fn_bf[2]))

    plt.plot(infile4[:,0], infile4[:,1], label='Time Step=0.25')
    plt.plot(infile5[:,0], infile5[:,1], label='Time Step=0.1')
    plt.plot(infile6[:,0], infile6[:,1], label='Time Step=0.05')
    plt.xlabel('Monte Carlo steps',fontsize=16)
    plt.ylabel(r'$E_L (\hbar \omega)$',fontsize=14)
    plt.show()
    return    





#variationalpha()
plottsteps()




