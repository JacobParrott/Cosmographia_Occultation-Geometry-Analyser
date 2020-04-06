
import numpy as np
import spiceypy as spice
import spiceypy.utils.support_types as stypes
import pandas as pd
from os import path
import matplotlib.pyplot as plt

#Import custom modules
import main
import atmosphere


#-----------------------------------------------------<VALUES TO EDIT REGULARLY>----------------------------------------
# If you only wish to analysis mutual [cross-link] occultation between MEX and TGO, then this is the only section that
# needs to be edited
start  = '2020 JAN 1'
stop   = '2020 JAN 3'
OCCSELECTION = 9 # Which occultation do you wish to see in Cosmographia? [optional]
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


# Inter the ingress epochs into a dataframe
occlist = np.ones((winsiz,3))
for i in range(winsiz):
    [ingress, egress] = spice.wnfetd(occwindow, i) # extract the begining and ends of the windows
    if i == 1 :
        ingresslist = ingress
    else:
        ingresslist = np.append(ingresslist, [ingress])


#form the dataframe
occs = pd.DataFrame(ingresslist, columns=['Time'])


#Form the profile geometry to be ported into Cosmographia (see README)
main.CosmographiaCatalogFormer(occs.Time[OCCSELECTION], sv)

#Calculate the residual doppler as the sum of the neutral and ionosphere
ray, dist = main.bending(occs.Time[OCCSELECTION], sv)
ionoresidual = atmosphere.iono(ray,dist,sv)
neutralresidual = atmosphere.neutral(ray,dist,sv)
residual = ionoresidual + neutralresidual

#Then integrate along residual to find total doppler


#Populate lists with geometric parameters for each epoch
lon = lat = dist = sza = angle = np.ones([winsiz,1])
for i in range(winsiz-1):

    #try to plot the location on a map with cartopy
    lon[i],lat[i], dist[i], nearestpoint = main.Location(occs.Time[i], sv)

    sza[i] = main.SolarZenithAngles(occs.Time[i],nearestpoint, sv)

    angle[i] = main.grazingangle(occs.Time[i], sv) 

    progess = (i/winsiz) *100
    print('%.2f' %progess)

#plot all of the tangent points onto the surface of mars
main.charter(lon, lat, start, stop,here)

# Add to dataframe
occs['Longitude'] = lon
occs['Latitude']  = lat
occs['Distance'] = dist
occs['SolarZenithAngle'] = sza
occs['GrazingAngle'] = angle










