# Common imports
import os
import numpy as np
import pandas as pd
from pandas import DataFrame
import matplotlib.pyplot as plt
#from scipy import *

CB91_Blue = '#2CBDFE'
CB91_Green = '#47DBCD'
CB91_Pink = '#F3A0F2'
CB91_Violet = '#661D98'
CB91_Amber = '#F5B14C'

color_list = [CB91_Blue, CB91_Pink, CB91_Green, CB91_Amber, CB91_Violet]
#plt.rcParams['axes.prop_cycle'] = plt.cycler(color=color_list)

def data_path(DATA_ID, dat_id):
    return os.path.join(DATA_ID, dat_id)

def variationalpha():
    """
    Plotting the results from gradient decent
    """
    #Filenames
    fn_noninteract=['03N10Dim3nice.txt', '07N10Dim3nice.txt', '03N100Dim3nice.txt']
    fn_interact=['N10Dim3new.txt', 'N50Dim3new.txt', 'N100Dim3new.txt']
    #Folders
    folder = ["Results/GDalpha/noninteract/bruteforce/analytic", "Results/GDalpha/interact/bruteforce/numeric"]

    infile = np.loadtxt(data_path(folder[1], fn_interact[0]))

    #plt.plot(infile[:,0], infile[:,1])
    #plt.plot(infile[:,0], infile[:,1], 'ro', markersize=3)
    #plt.xlabel(r'$\alpha$',fontsize=14)
    #plt.ylabel(r'$\langle E_L \rangle(\hbar \omega) $',fontsize=14)
    #plt.grid()
    #plt.show()

    # Generate data for the zoomed portion
    #X_detail=infile[5:14,0]
    #Y_detail = infile[5:14,1]
    #X_detail=infile[29:len(infile),0]
    #Y_detail = infile[29:len(infile),1]
    X_detail=infile[8:len(infile),0]
    Y_detail = infile[8:len(infile),1]

    # plot the main figure
    plt.plot(infile[:,0], infile[:,1])
    plt.plot(infile[:,0], infile[:,1], 'ro', markersize=3)
    plt.grid()
    plt.xlabel(r'$\alpha$',fontsize=14)
    plt.ylabel(r'$\langle E_L \rangle(\hbar \omega) $',fontsize=14)

    # location for the zoomed portion 
    sub_axes = plt.axes([.52, 0.57, 0.35, 0.25]) 
    #sub_axes = plt.axes([0.2, 0.57, 0.35, 0.25]) 
    #sub_axes = plt.axes([0.22, 0.59, 0.35, 0.25]) 

    # plot the zoomed portion
    sub_axes.plot(X_detail, Y_detail) 
    sub_axes.plot(X_detail, Y_detail, 'ro', markersize=3) 

    # insert the zoomed figure
    #plt.setp(sub_axes)
    plt.grid()

    plt.show()


    return


def linspacealpha():
    """
    Plotting the local energy as variation of alpha
    """
    #Dim=3, particles=10, noninteracting, 2^18 steps
    a_nonint=  [0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7]
    el_nonint= [16.98, 15.92, 15.38, 15.09, 15, 15.04, 15.26, 15.47, 15.92]

    plt.plot(a_nonint, el_nonint)
    plt.plot(a_nonint, el_nonint, 'ro', markersize=3)
    plt.xlabel(r'$\alpha$',fontsize=14)
    plt.ylabel(r'$\langle E_L \rangle(\hbar \omega) $',fontsize=14)
    plt.grid()
    plt.show()

    #Dim=3, particles=10, noninteracting, 2^18 steps
    a_is_nonint=  [0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7]
    el_is_nonint= [21.7391, 16.57, 15.203, 15.0384 ,15, 15.638, 16.387, 17.3116, 18.5878]

    plt.plot(a_is_nonint, el_is_nonint)
    plt.plot(a_is_nonint, el_is_nonint, 'ro', markersize=3)
    plt.xlabel(r'$\alpha$',fontsize=14)
    plt.ylabel(r'$\langle E_L \rangle(\hbar \omega) $',fontsize=14)
    plt.grid()
    plt.show()

    #Dim=3, particles=10, interacting, 2^16 steps
    a_int=  [0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7]
    el_int= [28.13, 25.76, 24.83, 24.55, 24.39, 24.54, 24.72, 25.15, 25.79]

    plt.plot(a_int, el_int)
    plt.plot(a_int, el_int, 'ro', markersize=3)
    plt.xlabel(r'$\alpha$',fontsize=14)
    plt.ylabel(r'$\langle E_L \rangle(\hbar \omega) $',fontsize=14)
    plt.grid()
    plt.show()

    return


def plottsteps():
    #Both plotted with alpha=0.4

    #Filenames
    fn_bf=['steplength005N10Dim3.txt', 'steplength01N10Dim3.txt', 'steplength025N10Dim3.txt', 'steplength05N10Dim3.txt', 'steplength1N10Dim3.txt']
    fn_is=['timestep0005N10Dim3.txt', 'timestep001N10Dim3.txt', 'timestep005N10Dim3.txt', 'timestep01N10Dim3.txt','timestep025N10Dim3.txt', 'timestep05N10Dim3.txt']

    #Folders
    folder= ["Results/steps/bruteforce/", "Results/steps/importancesampling/"]

    #infile = open(data_path(folder_noninteract[0], fn[4]),'r')

    infile = np.loadtxt(data_path(folder[0], fn_bf[0]))
    infile2 = np.loadtxt(data_path(folder[0], fn_bf[1]))
    infile3 = np.loadtxt(data_path(folder[0], fn_bf[2]))
    infile4 = np.loadtxt(data_path(folder[0], fn_bf[3]))
    infile5 = np.loadtxt(data_path(folder[0], fn_bf[4]))

    #plt.plot(infile[:,0], infile[:,1], label='Step Length=0.05')
    plt.plot(infile2[:,0], infile2[:,1], label='Step Length=0.1')
    plt.plot(infile3[:,0], infile3[:,1], label='Step Length=0.25')
    plt.plot(infile4[:,0], infile4[:,1], label='Step Length=0.5')
    plt.plot(infile5[:,0], infile5[:,1], label='Step Length=1')
    plt.xlabel('Monte Carlo cycles',fontsize=12)
    plt.ylabel(r'$\langle E_L \rangle (\hbar \omega)$',fontsize=14)
    plt.grid()
    plt.legend()
    plt.show()

    infile6 = np.loadtxt(data_path(folder[1], fn_is[0]))
    infile7 = np.loadtxt(data_path(folder[1], fn_is[1]))
    infile8 = np.loadtxt(data_path(folder[1], fn_is[2]))
    infile9 = np.loadtxt(data_path(folder[1], fn_is[3]))
    infile10 = np.loadtxt(data_path(folder[1], fn_is[4]))
    infile11 = np.loadtxt(data_path(folder[1], fn_is[5]))

    #plt.plot(infile6[:,0], infile6[:,1], label='Time Step=0.005')
    plt.plot(infile7[:,0], infile7[:,1], label='Time Step=0.01')
    plt.plot(infile8[:,0], infile8[:,1], label='Time Step=0.05')
    plt.plot(infile9[:,0], infile9[:,1], label='Time Step=0.1')
    plt.plot(infile10[:,0], infile10[:,1], label='Time Step=0.25')
    #plt.plot(infile11[:,0], infile11[:,1], label='Time Step=0.5')

    plt.xlabel('Monte Carlo cycles',fontsize=12)
    plt.ylabel(r'$\langle E_L \rangle(\hbar \omega) $',fontsize=14)
    plt.grid()
    plt.legend(loc=(0.67, 0.67))
    plt.show()
    
    return



variationalpha()
#plottsteps()
#linspacealpha()



