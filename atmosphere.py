#I/O Input is the positions of the two space craft 
# output is firstly the bending angle, then expected doppler
# for all atmospheric effects look at https://academic.oup.com/astrogeo/article/44/4/4.6/195844#1924741

import numpy as np
import spiceypy as spice
import math
from scipy import constants
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.optimize import curve_fit 

# from mpmath import *
# mp.dps = 30 ; mp.pretty = True 

#Equation defining the Chapman layer model to describe refrective indecies through the martian ionosphere. 
#Equation 1 & constants from https://descanso.jpl.nasa.gov/propagation/mars/MarsPub_sec2.pdf.
#THIS FUNCTION IS NOT CORRECTED FOR SZA YET
def iono(ray,dist):
    alt =  1000*ray #CONVERT TO M
    #sza = ray[3,:]
    #sza = np.zeros([1,dist])
    Ne1 = np.zeros([1,dist]) 
    Ne2 = np.zeros([1,dist])#MATRIX
    N1 = np.zeros([1,dist])
    N2 = np.zeros([1,dist])
    
    H = 13000 #scale height for ionosphere
    
    h = 135000 # Peak N's altitude [used to be 130000km]
    N0 = 1.2e11 #Peak electron density (FOR REALISTIC PLOTS 1e11 )
    
    # need to treat the SZA vector so that it doesnt exceed. (THIS WILL ONLY WORK DURING DAY OCCULTATIONS)

    sza = 0
    for i in range(dist):
        #densityadjustment = N0 * (np.cos(sza[i]) ** k ) # SZA adjustment
        NHP = (alt[i] - h) / H #Normalised Height Parameter
        term1 =(1/ np.cos(sza)) * np.exp(-NHP)
        exponetial = (1 - NHP - (term1))/2
        Ne1[0,i] =  N0 * np.exp(exponetial)
        N1[0,i] = (-40.2592397908339) * (Ne1[0,i] / (437.1e6 ** 2))


    #Ao. et al has another bump slightly lower, these needs to be replicated
    _H = 8000 #scale height for ionosphere
    _h = 100000 # Peak N's altitude
    _N0 = 3e10 #Peak electron density (FOR REALISTIC PLOTS 1e11 )

    for i in range(dist):
        #densityadjustment = N0 * (np.cos(sza[i]) ** k ) # SZA adjustment
        NHP = (alt[i] - _h) / _H #Normalised Height Parameter
        term1 =(1/ np.cos(sza)) * np.exp(-NHP)
        exponetial = (1 - NHP - (term1))/2
        Ne2[0,i] =  _N0 * np.exp(exponetial)
        N2[0,i] = (-40.2592397908339) * (Ne2[0,i] / (437.1e6 ** 2))
    
    Ne = Ne1+Ne2
    N =N1+N2
    return N, Ne #^this might be wrong (^-3) only done to make a graph look right





def neutral(ray,dist):#NEEDS A LIMIT, IF ABOVE 100KM, CUT TO 0
    alt = 1000* ray
    Nn = np.ones([1,dist])
    N = np.ones([1,dist])
    pres = np.zeros([1,dist])
    denom = np.zeros([1,dist])
    #Martian height and pressure have different models depending on alt
    for i in range(dist):
            temp = (-23.4 - (9.98e-4 * alt[i])) +273.1 # convert to kelvin
            if temp<0:
                temp=0 #you cannot have negative kelvin
                Nn[0,i] = np.inf
                continue
            pres[0,i] = 610 * np.exp(-9e-5 * alt[i])
            denom[0,i] = (constants.Boltzmann * temp) #increase the temperature, density reduces

            Nn[0,i] = pres[0,i]/denom[0,i]

    # if we go negative then we should set to 0, this model is far from perfect
    Nn[Nn == np.inf] = 0
    N = Nn * 1.804e-29 #(1.804e-29)
    #Nn[0,101:113] = 0
    return N,Nn


def new_neutral(ray, dist):
    H = 11
    Max_N = 3.9e-6
    alt =  ray
    N = np.ones([1,dist])
    for i in range(dist):
        N[0,i] = Max_N * np.exp(-alt[i]/H)
    return N

# def exponential(x, a, b):
#     return a*np.exp(b*x)


# CAN PLOT THE SHAPE BELOW, SHOULD LOOK LIKE WHAT BRUNO HAD
#lets produce a lil smoothed analytic of the nEurtral
height = np.arange(0,300,0.1) 
ray = np.zeros(len(height))
ray[:] = height
# lets model to 300 km 
NeutralIndex = new_neutral(ray, len(height))
IonoIndex, ElectronDensity = iono(ray, len(height))
net_N = ((NeutralIndex) + IonoIndex) *1e6
ElectronIndex, ElectronDesity = iono(ray, len(height))


# #we need a derivate of the iono density

smoothresults = InterpolatedUnivariateSpline(height, ElectronDesity)
prime = smoothresults.derivative(n=1)
doubleprime = smoothresults.derivative(n=2)
finerx = np.linspace(0,300,1000)
firstderivative= prime(finerx)
secondderivative= doubleprime(finerx)

        # fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2)

        # ax1.plot(NeutralDesity[0,:],height, 'b')
        # ax1.set_title('Neutral Density', loc='left')
        # ax1.set_ylabel('Altitude (km)')
        # ax1.set_xlabel('Density ($m^{-3}$)')
        # ax1.set_ylim(0,300)

        # ax2.plot(np.log10(NeutralDesity[0,:]),height,'b')
        # ax2.set_title('Neutral Density', loc='left')
        # ax2.set_ylabel('Altitude (km)')
        # ax2.set_xlabel(' LOG Density ($m^{-3}$)')
        # ax2.set_ylim(0,300)

        # ax3.plot(ElectronDesity[0,:],height,'m')
        # ax3.set_title('Electron Density', loc='left')
        # ax3.set_ylabel('Altitude (km)')
        # ax3.set_xlabel('Density ($m^{-3}$)')
        # ax3.set_ylim(0,300)


        # ax4.plot(np.log10(ElectronDesity[0,:]),height,'m')
        # ax4.set_title('Electron Density', loc='left')
        # ax4.set_ylabel('Altitude (km)')
        # ax4.set_xlabel(' LOG Density ($m^{-3}$)')
        # ax4.set_ylim(0,300)

        # plt.show()


plt.subplot(121)
plt.plot(net_N[0,:],height, '-.')
plt.title('Ionospheric + Atmospheric Refractivty', loc='left')
plt.ylabel('Altitude (km)')
plt.xlabel('Refractivity (dimentionless)')

plt.subplot(122)
plt.plot(-firstderivative,finerx)
plt.title(' Gradient in Electron Density', loc='left')
plt.ylabel('Altitude (km)')
plt.xlabel('Delta Density ($m^{-3}$)')

# #ax[2].plot(-firstderivative,finerx)



plt.show()



# # plt.rcParams['text.usetex']=True
# # plt.rcParams['text.latex.unicode']=True

# plt.plot( N[0], height)
# plt.title('A Plot to Show Refractive Index Through the Martian Neutral Atmosphere')
# plt.xlabel(r'$\Delta$N, Refractive Index' )
# plt.ylabel('Altitude, using a spherical model for Mars (km)')
# plt.show()



# ynew = np.ones([1,300])
# x = height
# y = N[0]
# a = np.polyfit(x,y,10)
# for i in height:
#     ynew[0,i] = np.polyval(a,x[i])







# fig,ax1 = plt.subplots()
# ax1.plot(N[0], height, linewidth = 5, color = 'y')
# #lets blur the dicontinuity
# #smoothresults = UnivariateSpline(height, N,s=0,k=5)
# ax2 = ax1.twinx()
# ax2.plot( ynew[0,:],height, ':')
# plt.show()


