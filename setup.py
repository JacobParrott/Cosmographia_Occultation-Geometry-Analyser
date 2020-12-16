

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
import math

#Import custom modules

import main
import atmosphere
import swiftmain #a faster version for the costly doppler profile calculations
import FindDopplerMain




#-----------------------------------------------------<VALUES TO EDIT REGULARLY>----------------------------------------
# If you only wish to analysis mutual [cross-link] occultation between MEX and TGO, then this is the only section that
# needs to be edited
start  = '2020 NOV 1'
stop   = '2020 NOV 5'
OCCSELECTION = 2 # Which occultation do you wish to see in Cosmographia? [optional]
here = path.abspath(path.dirname(__file__))

PathtoMetaKernel1 = here + '/TGO/mk/em16_plan.tm'
PathtoMetaKernel2 = here + '/MEX/krns/mk/MEX_OPS.tm'

print(PathtoMetaKernel1)
print(PathtoMetaKernel2)

#experiment with the old and new kernels (hopfully they are the exact same)
# PathtoMetaKernel1 = here + '/TGOold/krns/mk/em16_plan.tm'
# PathtoMetaKernel2 = here + '/MEXold/krns/mk/MEX_OPS.tm'

#PathtoMetaKernel3 = here + '/MarsRec_Odyssey/MRO_MO_mk.tm'
#-----------------------------------------------------------------------------------------------------------------------


spice.furnsh(PathtoMetaKernel1)
spice.furnsh(PathtoMetaKernel2)
#spice.furnsh(PathtoMetaKernel3)

sv = main.SpiceVariables()

#~~~~~~~EXPERIMENTZONE#######################
from scipy import constants
import pickle
geometricdopplershift = np.zeros(600)
for time in range(599,0,-1):
    sc2scstates = spice.spkezr(sv.target, (657605289.1825405 - time), sv.fframe, 'LT+S', sv.obs)
    velocityvector = sc2scstates[0][3:6]
    velocityvector = velocityvector[0:3]
    positionalvector =  sc2scstates[0][0:3]
    positionalvector = positionalvector[0:3]
    velocityangle = spice.vsep( positionalvector, velocityvector) #rads
    relativevelocity = np.linalg.norm(velocityvector) * np.cos(velocityangle) 
    
    geometricdopplershift[time] = -(relativevelocity/constants.c) * 437.1e9 #conversion from km to m

pd.DataFrame(geometricdopplershift).to_csv("geometricdopplershift.csv")






# Setting Variables
ingresslist = np.array([1.0] ,dtype = float)
egresslist = np.array([1.0] ,dtype = float)
etbeg = spice.str2et(start)
etend = spice.str2et(stop)


# Form a windows that gfoclt can populate
window = stypes.SPICEDOUBLE_CELL(2)
spice.wninsd(etbeg, etend,window)
occwindow = stypes.SPICEDOUBLE_CELL(sv.MAXILV)

#find occultation windows between the dates listed above [ most comp cost in this function]
spice.gfoclt('ANY',sv.front,sv.fshape,sv.fframe, sv.target, sv.bshape, '' , sv.abcorr, sv.obs, sv.stepsz, window, occwindow)

winsiz = spice.wncard( occwindow )# Find cardinality (number of windows)

#initialize lists to form dataframe
lon , lat, dist, sza, angle = ( np.ones(winsiz-1) for i in range(5))


# Enter the ingress epochs into a dataframe
occlist = np.ones((winsiz,3))
for i in range(winsiz):
    [ingress, egress] = spice.wnfetd(occwindow, i) # extract the begining and ends of the windows
    if i == 1 :
        ingresslist = ingress
        egresslist = egress
    else:
        ingresslist = np.append(ingresslist, [ingress])
        egresslist = np.append(egresslist, [egress])

#MEX,TGO = FindDopplerMain.ephemerides(636391881,1)
#form the dataframe
occs = pd.DataFrame(ingresslist, columns=['Time'])
occs['Ingress'] = ingresslist
occs['Egress'] = egresslist


#Form the profile geometry to be ported into Cosmographia (see README)


# Produce geometry in wavelenghts to investigate electric distance

date_time = [""]*(winsiz-1)
ingress_date_time = [""]*(winsiz-1)
localtime = [""]*(winsiz-1)
lon = np.ones([winsiz-1,1])
veldopp = np.ones([winsiz-1,1])
lat =  np.ones([winsiz-1,1])
dist=np.ones([winsiz-1,1])
sza= np.ones([winsiz-1,1])
clearanceangle = np.ones([winsiz-1,1])
Rx = np.ones([winsiz-1,1])
MexNadirTGOAngle = np.ones([winsiz-1,1])
#grazingangle = np.ones([winsiz-1,1]) # HERE

for i in tqdm(range(winsiz-1)):

    #try to plot the location on a map with cartopy
    lon[i],lat[i], dist[i], nearestpoint, alt, speed, date_time[i],ingress_date_time[i], veldopp[i], MexNadirTGOAngle[i] = main.Location(occs.Time[i],occs.Ingress[i], sv, 0)# FUNCTION NEEDS NEW NAME
    sza[i] = main.SolarZenithAngles(occs.Time[i],nearestpoint, sv)
    clearanceangle[i] = main.earlyclearance(occs.Time[i],sv)
    hour, minute,sec,_,_ = spice.et2lst(occs.Time[i],499,(lon[i]* (math.pi/180)), 'PLANETOCENTRIC') #499 = spice code for Mars
    localtimevect = [hour, minute, sec]
    localtime[i]  = "%d:%d:%d" %(localtimevect[0],localtimevect[1],localtimevect[2])
    #angle[i] = main.grazingangle(occs.Time[i], sv) # this process is probably slow because this function is forming the entire profile.
    Rx[i] = main.expectedpower(occs.Time[i],sv)
    

    # progess = (i/(winsiz-1)) *100
    # print('%.2f' %progess)

#plot all of the tangent points onto the surface of mars
#main.charter(lon, lat, start, stop,here)

# Add to dataframe
occs['Julian'] = date_time
#occs['gress'] = inress_date_time
occs['Longitude (°E)'] = lon[:,0]
occs['Latitude(°N)']  = lat[:,0]
occs['Distance (km)'] = dist[:,0]
occs['VelocityDoppler'] = veldopp[:,0]
occs['SolarZenithAngle(°)'] = sza[:,0]
occs['LocalTime'] = localtime
occs['6MinPriorClearance (°)'] = clearanceangle[:,0]
occs['Rx (dBm)'] = Rx[:,0]
occs['Mex Nadir to TGO Angle'] = MexNadirTGOAngle[:,0]
#occs['GrazingAngle'] = angle[:] #currently excluded for speed reasons.


occs = occs[(occs["Distance (km)"] <10000)]# ADD FURTHER FILTERING HERE

occs.to_excel("OccultationSelection.xlsx")#Will not save to exel file if it is open in another window

CosmoTime = ingresslist[OCCSELECTION]
main.CosmographiaCatalogFormer(658560507, sv)


print('stop here')
#Example filtering function: nba[    (nba["_iscopy"] == 0) &    (nba["pts"] > 100) &    (nba["opp_pts"] > 100) &    (nba["team_id"] == "BLB")]



#THIS IS THE INTENSIVE REFRACTIVE BENDING SECTION
####################################################################################################
#ALL OF THIS MIGHT BE PUSHED TO A NEW PROGRAMME AS ITS PURPOSE HAS DIVERGED FROM THE ORIGINAL PURPOSE (STILL USEFUL THOUGH)
#ray, dist, totalperiods, remainingdistance = main.producegeometrylamda(occs.Time[OCCSELECTION], sv, 5)
# Ssize = 1
# S = np.zeros(Ssize)

# for time in range(Ssize):
#     initialangle, MEX,TGO, xyzpoints= main.producegeometrymeter(631255812.2432083,sv, time)
#     Bending, ElectricDistance = main.flatbending(xyzpoints, initialangle, sv, MEX, TGO)
#     S[time] = ElectricDistance
 
# np.savetxt('ElectricDistance.csv', S, delimiter = ',')
# print("this worked")

# with open('ElectricDistance.csv') as Scsv:
#     Snew = list(csv.reader(Scsv)) #record all the electric distance values in km

# ####################################################################################################













