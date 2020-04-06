#I/O Input is the positions of the two space craft 
# output is firstly the bending angle, then expected doppler

import numpy as np
import spiceypy as spice
import math
from scipy import constants

#Equation defining the Chapman layer model to describe refrective indecies through the martian ionosphere. 
#Equation 1 & constants from https://descanso.jpl.nasa.gov/propagation/mars/MarsPub_sec2.pdf.
#THIS FUNCTION IS NOT CORRECTED FOR SZA YET
def iono(ray,dist,sv):
    alt = 1000 * ray[2,:]
    #sza = ray[3,:]
    #sza = np.zeros([1,dist])
    Ne = np.ones([1,dist])
    N = np.ones([1,dist])
    H = 13000 #scale height for ionosphere
    k = 0.57 #magic number
    h = 125000 # Peak N's altitude
    N0 = 2e5 #Peak electron density
    
    # need to treat the SZA vector so that it doesnt exceed. (THIS WILL ONLY WORK DURING DAY OCCULTATIONS)

    sza = 0
    for i in range(dist):
        #densityadjustment = N0 * (np.cos(sza[i]) ** k ) # SZA adjustment
        NHP = (alt[i] - h) / H #Normalised Height Parameter
        if i==10000:
            print('stop here')
        thing =(1/ np.cos(sza)) * np.exp(-NHP)
        exponetial = 1 - NHP - (thing)
        Ne[0,i] = 1e6 * N0 * np.exp( exponetial)
        N[0,i] = (-40.31) * (Ne[0,i] / (437.1e6 ** 2)) #convert electron density to refractive index {-k * Ne/f^2}
    
    return N


def neutral(ray,dist,sv):
    alt = 1000 * ray[2,:]
    Nn = np.ones([1,dist])
    N = np.ones([1,dist])
    #Martian height and pressure have different models depending on alt
    for i in range(dist):
        if alt[i] >= 7000:
            temp = -23.4 - (0.00222 * alt[i])
        else:
            temp = -31 -(9.98e-4 * alt[i])

        pres = 699 * np.exp(-9e-5 * alt[i])
        density = pres / (0.1921 *(temp +273.1))
        Nn[0,i] = pres/(constants.Boltzmann * (temp +273.1))
        N[0,i] = Nn[0,i] * 1.804e-29

    return N    


