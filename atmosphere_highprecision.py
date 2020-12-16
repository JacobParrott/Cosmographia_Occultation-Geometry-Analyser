#I/O Input is the positions of the two space craft 
# output is firstly the bending angle, then expected doppler

import numpy as np
import spiceypy as spice
import math
from scipy import constants
from mpmath import *
mp.dps = 60 ; mp.pretty = True 
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
from scipy.optimize import curve_fit 

#Equation defining the Chapman layer model to describe refrective indecies through the martian ionosphere. 
#Equation 1 & constants from https://descanso.jpl.nasa.gov/propagation/mars/MarsPub_sec2.pdf.
#THIS FUNCTION IS NOT CORRECTED FOR SZA YET
def iono(ray,dist):
    alt = 1000 * ray #ELEECTRIC DISTANCE CALCS NEED 1000S OF THESE, BUT UR BENDING IS IN KM
    #sza = ray[3,:]
    #sza = np.zeros([1,dist])
    Ne = matrix(1,dist)
    N = matrix(1,dist)
    H = mpf(15000) #scale height for ionosphere
    k = mpf(0.57) #magic number
    h = mpf(130000) # Peak N's altitude
    N0 = mpf(2.000000000000001e11) #Peak electron density ALTERATION HERE (Mars =2e11)
    
    # need to treat the SZA vector so that it doesnt exceed. (THIS WILL ONLY WORK DURING DAY OCCULTATIONS)

    #**********************************************
    #GET RID OF ALL DECIMAL PACK USAGE AND REPLACE WITH MPMATH. COLLAPSE ALL SUBLINES INTO ONE TO BE CONSISE
    #**********************************************
    sza = 0
    for i in range(dist):
        ALT = mpf(alt[i])
        #densityadjustment = N0 * (np.cos(sza[i]) ** k ) # SZA adjustment
        NHP = ((ALT) - h) / H #Normalised Height Parameter
        if i==10000:
            print('stop here')
        thing =(1/ np.cos(sza)) * (math.e ** (-NHP))
        exponetial = 1 - NHP - (thing)
        otherthing = (math.e ** exponetial)
        densitystore = N0 * otherthing
        Ne[0,i] = densitystore
        N[0,i] =  -40.2592397908339 * (Ne[0,i] / 437.1e6 ** 2) #convert electron density to refractive index {-k * Ne/f^2}
    
 
    return N





def neutral(ray,dist):
    
    alt = 1000 * ray
    
    Nn = matrix(1,dist)
    N = matrix(1,dist)
   
    # = np.zeros([1,dist])
    #Martian height and pressure have different models depending on alt
    for i in range(dist):

        ALT = mpf(alt[i])
        temp = -23.4 -(9.98e-4 * ALT) + 273.1
        if temp<0:
            temp=0 #you cannot have negative kelvin
            Nn[0,i] = 0
            continue
        pres = mpf(610) * exp(-9e-5 * ALT) 
        Nn[0,i] = pres/(constants.Boltzmann * temp)
        N[0,i] = Nn[0,i] * 1.804e-29 
        #N[0,i] = temp


    return N 



