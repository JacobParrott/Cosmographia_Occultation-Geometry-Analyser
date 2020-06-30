

import numpy as np
import spiceypy as spice
import spiceypy.utils.support_types as stypes
import pandas as pd
from os import *
import matplotlib.pyplot as plt
import time
from tqdm import tqdm
import csv
from multiprocessing import Pool

#Import custom modules

import main
import atmosphere
import swiftmain #a faster version for the costly doppler profile calculations
import FindDopplerMain




#-----------------------------------------------------<VALUES TO EDIT REGULARLY>----------------------------------------
# If you only wish to analysis mutual [cross-link] occultation between MEX and TGO, then this is the only section that
# needs to be edited
start  = '2020 JAN 1'
stop   = '2020 JAN 3'
OCCSELECTION = 14 # Which occultation do you wish to see in Cosmographia? [optional]
here = path.abspath(path.dirname(__file__))
PathtoMetaKernel1 = here + '/TGO/krns/mk/em16_plan.tm'
PathtoMetaKernel2 = here + '/MEX/krns/mk/MEX_OPS.tm'
#-----------------------------------------------------------------------------------------------------------------------

spice.furnsh(PathtoMetaKernel1)
spice.furnsh(PathtoMetaKernel2)

sv = main.SpiceVariables()

# Setting Variables
ingresslist = np.array([1.0] ,dtype = float)
etbeg = spice.str2et(start)
etend = spice.str2et(stop)


# Form a windows that gfoclt can populate
window = stypes.SPICEDOUBLE_CELL(2)
spice.wninsd(etbeg, etend,window)
occwindow = stypes.SPICEDOUBLE_CELL(sv.MAXILV)

#find occultation windows between the dates listed above [ most comp cost in this function]
spice.gfoclt('ANY',sv.front,sv.fshape,sv.fframe, sv.target, sv.bshape, 'J2000' , sv.abcorr, sv.obs, sv.stepsz, window, occwindow)

winsiz = spice.wncard( occwindow )# Find cardinality (number of windows)

#initialize lists to form dataframe
lon , lat, dist, sza, angle = ( np.ones(winsiz-1) for i in range(5))


# Enter the ingress epochs into a dataframe
occlist = np.ones((winsiz,3))
for i in range(winsiz):
    [ingress, egress] = spice.wnfetd(occwindow, i) # extract the begining and ends of the windows
    if i == 1 :
        ingresslist = ingress
    else:
        ingresslist = np.append(ingresslist, [ingress])

#MEX,TGO = FindDopplerMain.ephemerides(636391881,1)
#form the dataframe
occs = pd.DataFrame(ingresslist, columns=['Time'])
print(occs.Time[OCCSELECTION])

#Form the profile geometry to be ported into Cosmographia (see README)
#main.CosmographiaCatalogFormer(occs.Time[OCCSELECTION], sv)

#Calculate the residual doppler as the sum of the neutral and ionosphere

# Produce geometry in wavelenghts to investigate electric distance

#ray, dist, totalperiods, remainingdistance = main.producegeometrylamda(occs.Time[OCCSELECTION], sv, 5)
Ssize = 1
S = np.zeros(Ssize)

for time in range(Ssize):
    initialangle, MEX,TGO, xyzpoints= main.producegeometrymeter(631255812.2432083,sv, time)
    Bending, ElectricDistance = main.flatbending(xyzpoints, initialangle, sv, MEX, TGO)
    S[time] = ElectricDistance
 
np.savetxt('ElectricDistance.csv', S, delimiter = ',')
print("this worked")

with open('ElectricDistance.csv') as Scsv:
    Snew = list(csv.reader(Scsv)) #record all the electric distance values in km




#####################################################################################
toc = time.clock()
print(toc-tic)
ionoresidual = atmosphere.iono(ray[2,:],totalperiods)
neutralresidual = atmosphere.neutral(ray[2,:],totalperiods)
residual = 1 + (ionoresidual + neutralresidual)

plt.plot(range(totalperiods),residual[0,:])
plt.title("Refractive Index through Propergation of Ray")
plt.xlabel("MEX->TGO distance (km)")
plt.ylabel("Refractive Index")
plt.show()

[electricdistance, geometric, dopplershift] = main.doppler(residual, totalperiods, dist, remainingdistance)
miss = electricdistance - dist + remainingdistance # a possitive number as electric distance is greater that geometric due to iono
howwrongurcodeis = geometric - dist # is this purly due to rounding of that 9945 (each 1 is ~685 m)
# print("Delta between electric distance and geometric distance is", "{:.8f}".format(miss * 1000 ), "m")

#Populate lists with geometric parameters for each epoch
lon = np.ones([winsiz,1])
lat =  np.ones([winsiz,1])
dist=np.ones([winsiz,1])
sza= np.ones([winsiz,1])
angle = np.ones([winsiz,1]) # HERE

for i in range(winsiz-1):

    #try to plot the location on a map with cartopy
    lon[i],lat[i], dist[i], nearestpoint, alt = main.Location(occs.Time[i], sv, 0)

    sza[i] = main.SolarZenithAngles(occs.Time[i],nearestpoint, sv)

    angle[i] = main.grazingangle(occs.Time[i], sv) 

    progess = (i/winsiz) *100
    print('%.2f' %progess)

#plot all of the tangent points onto the surface of mars
main.charter(lon, lat, start, stop,here)

# Add to dataframe
occs['Longitude'] = lon[0,:]
occs['Latitude']  = lat[0,:]
occs['Distance'] = dist[0,:]
occs['SolarZenithAngle'] = sza[0,:]
occs['GrazingAngle'] = angle[0,:]










