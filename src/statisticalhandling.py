"""
Code on blocking from the article of Marius Jonsson
"""
# Common imports
import os


def data_path(DATA_ID, dat_id):
    return os.path.join(DATA_ID, dat_id)

from numpy import log2, zeros, mean, var, sum, loadtxt, arange, array, cumsum, dot, transpose, diagonal, sqrt
from numpy.linalg import inv

def block(x):
    # preliminaries
    n = len(x)
    print(n)
    d = int(log2(n))
    s, gamma = zeros(d), zeros(d)
    mu = mean(x)

    # estimate the auto-covariance and variances 
    # for each blocking transformation
    for i in arange(0,d):
        n = len(x)
        # estimate autocovariance of x
        gamma[i] = (n)**(-1)*sum( (x[0:(n-1)]-mu)*(x[1:n]-mu) )
        # estimate variance of x
        s[i] = var(x)
        # perform blocking transformation
        #print(x[0::2])
        #print(x[1::2])
        x = 0.5*(x[0::2] + x[1::2])
   
    # generate the test observator M_k from the theorem
    M = (cumsum( ((gamma/s)**2*2**arange(1,d+1)[::-1])[::-1] )  )[::-1]

    # we need a list of magic numbers
    q =array([6.634897,9.210340, 11.344867, 13.276704, 15.086272, 16.811894, 18.475307, 20.090235, 21.665994, 23.209251, 24.724970, 26.216967, 27.688250, 29.141238, 30.577914, 31.999927, 33.408664, 34.805306, 36.190869, 37.566235, 38.932173, 40.289360, 41.638398, 42.979820, 44.314105, 45.641683, 46.962942, 48.278236, 49.587884, 50.892181])

    # use magic to determine when we should have stopped blocking
    for k in arange(0,d):
        if(M[k] < q[k]):
            break
    if (k >= d-1):
        print("Warning: Use more data")
    return mu, s[k]/2**(d-k)

#Defining the filenames
fn=   ["N=1Dim=1", "N=10Dim=1", "N=100Dim=1", "N=500Dim=1"
        ,"N=1Dim=2", "N=10Dim=2", "N=100Dim=2", "N=500Dim=2"
        ,"N=1Dim=3", "N=10Dim=3", "N=100Dim=3", "N=500Dim=3"]


#Folders noninteracting
folder_bruteforce = ["Results/bruteforce/analytic", "Results/bruteforce/numeric"]
folder_importance = ["Results/importancesampling/analytic", "Results/importancesampling/numeric"]

#Filenames interacting
fn_interact=   ["N=10Dim=3_interact_new", "N=50Dim=3_interact", "N=100Dim=3_interact", "N=10Dim=3_interact_220steps"]
folder_interact = ["Results/bruteforce/numeric/interact/"]

infile = open(data_path(folder_bruteforce[1], fn[3]),'r')

x = loadtxt(infile)
(mean, var) = block(x) 
std = sqrt(var)
import pandas as pd
from pandas import DataFrame
data ={'Mean':[mean], 'STDev':[std]}
frame = pd.DataFrame(data,index=['Values'])
print(frame)