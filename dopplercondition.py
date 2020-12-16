#dopplercondition
import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy import constants
from scipy import signal
import FindDopplerMain
import main
from os import path
import spiceypy as spice


#need to call it once to see if the distance has regular jumps in it
#-----------------------------------------------------<VALUES TO EDIT REGULARLY>----------------------------------------
# If you only wish to analysis mutual [cross-link] occultation between MEX and TGO, then this is the only section that
# needs to be edited
start  = '2020 JAN 1'
stop   = '2020 JAN 3'
OCCSELECTION = 14 # Which occultation do you wish to see in Cosmographia? [optional]
here = path.abspath(path.dirname(__file__))
PathtoMetaKernel1 = here + '/TGO/mk/em16_plan.tm'
PathtoMetaKernel2 = here + '/MEX/mk/MEX_OPS.tm'
#-----------------------------------------------------------------------------------------------------------------------

spice.furnsh(PathtoMetaKernel1)
spice.furnsh(PathtoMetaKernel2)

sv = main.SpiceVariables()
# initialangle, MEX,TGO, xyzpoints= main.producegeometrymeter(636491202.20059,sv, 3)
# Bending, ElectricDistance = main.flatbending(xyzpoints, initialangle, sv, MEX, TGO)

# S = []
# GeoDistance = []
# #SORTED, NOW YOU NEED TO CALL THIS FUNCTION 60 TIMES 
# # AS EPHEMERIDES WILL ONLY GIVE A STRIAGHT LINE, THEN SEE IF 
# # SPICE WILL LET YOU SUBSAMPLE, MAYBE THAT WONT HAVE TEETH
# for i in np.arange(60,0,-1):
#    ElectricDistance, Geometric  = FindDopplerMain.thecall(i)
# # #S2, GeoDistance2  = FindDopplerMain.thecall(2)
#    S = np.append( S,ElectricDistance ) ; GeoDistance = np.append(GeoDistance, Geometric)

#CONTROL TEST
#S, GeoDistance  = FindDopplerMain.thecall(3)


results = []

with open('ElectricDistance.csv') as csvDataFile:
   csvReader = csv.reader(csvDataFile)
   for row in csvReader:
      results.append(row)

Results = np.array(results)

[ylen, xlen] = Results.shape



s = Results[0,:]
geodistance = Results[1,:]

tangentalt = Results[2,:]


S= list(map(float, s))
GeoDistance= list(map(float, geodistance))
TangentAlt= list(map(float, tangentalt))

S = np.array(S)
GeoDistance = np.array(GeoDistance)
TangentAlt = (np.array(TangentAlt))/1000
#S = np.flip(np.array(S))  #muliprocessing produces the dataset backwards
#GeoDistance = np.flip(np.array(GeoDistance))
#TangentAlt = np.flip(np.array(TangentAlt))

#Results = np.array(results)) 
#results = np.flip(results)

#create a mask filter [1 means remove]

      # maskfilter = np.ones((xlen,ylen))
      # for event in range(xlen):
      #    if Results[event, 2] > 0:
      #       maskfilter[event, :] = 0
      # Results = np.ma.masked_array(Results, maskfilter)

print('results = ', results)
print('S = ', S)
print('GeoDistance = ', GeoDistance)
height = (np.size(results,0))
width = (np.size(results,1))

print('width = ',width)
print('height = ',height)
datasetsize = height * width #should give 50


#convert to float to perfrom calcs on
# a = np.zeros(datasetsize) 
# p=0
# for j in range(height):#move through the seconds count
#    for i in range(width): # move through the Hz count
#       value = S[j][i]
#       a[p] = float(value)
#       p = p+1

#print(a[1:width])

minA = np.amin(S)# rounding is ok if u are not iterating
mingeo = np.amin(GeoDistance)
#print(a)
print(minA)
#print('array = ',a)
ConditionedS = S - minA#convert to m from km
GeoDistance= GeoDistance-mingeo 


#for i in range(datasetsize):
   # ConditionedS[i] =ConditionedS[i] % 1 #remove the step function, BUT WHY DO WE HAVE A STEP FUNCTION?? THE KEY

print('MODarray = ',ConditionedS)
# this needs to be dynamic to the size of dataset. THIS IS MESSY

#ConditionedS = ConditionedS[beginningcut:] #removing the beginning of the dataset
# x = np.arange(0.0,(height - (beginningcut/10)),0.1) #4.4
# dx = np.arange(0.0,(height - ((beginningcut+1)/10)),0.1) #4.3

x = np.arange(0.0,len(S),1) #4.4
dx = np.arange(0.0,(len(S))-1,1) #4.3

# x = np.arange(len(S),0.0,-1) #4.4
# dx = np.arange((len(S))-1,0.0,-1) #4.3


#extract the dervative of S (x increments uniform so dx =1)
d =  np.diff(S)
dgeo = np.diff(GeoDistance)
dgeo = dgeo * (437.1e9/constants.c)
d = d * (437.1e9/constants.c) #should be in MHz, but we convert m to km so * 1000

#NEED TO REMOVE THE SUB-SECOND ALLIASING, IMPLEMENT A SMOOTHING FUNCTION 
for i in range(5, len(d)):
   rollingavg = (sum(d[i-4:i-1])/3)
   print((d[i-4:i-1]))
   #rollingavggeo = (sum(dgeo[i-4:i-1])/3) 
   if (d[i] < (rollingavg - 30)): 
      d[i] = rollingavg 
   elif (d[i] > (rollingavg + 30)):
      d[i] = rollingavg 

#applying a butterworth filter to avoid shot noise in rotation error
# period = len(d)
# maxfreq = 100
# nyquist = 0.5
# num,denom = signal.butter(2, 0.05, btype = 'low', analog = False)
# d = signal.filtfilt(num, denom, d)


#make a moving average whilst my sampling frequency is too low
# for i in range(np.size(dx)):
#     if i> np.size(dx)-3:
#         break
#     d[i] = d[i] + d[i+1]+ d[i+2]

residual = d-dgeo

greycolour = '#E3D9D7'


fig, ax = plt.subplots(2,2)
ax[0][0].set_title('Normalised Electric Distance and Geometric Distance (km)')
ax[0][0].plot(x,ConditionedS, 'r-')
ax[0][0].plot(x,GeoDistance, 'g*' )
ax[0][0].set_xlabel('Time (s) [end = occultation epoch]')
ax[0][0].set_ylabel('S (km) [:]', color = 'red')
#ax02 = ax[0][0].twinx()
#ax02.plot(x,TangentAlt, greycolour,linestyle =  ':' )

ax[0][1].plot(dx,d, 'g')
ax[0][1].plot(dx,dgeo,'r')
ax[0][1].set_title('Total Doppler Shift [residual + geometric](Hz)')
#ax12 = ax[0][1].twinx()
#ax12.plot(x,TangentAlt, greycolour, linestyle = ':' )
ax[0][1].grid(True)

ax[1][0].set_title('LEFT BLANK')
#ax[1][0].plot(dx,dgeo,'r')
ax[1][0].grid(True)
#ax32 = ax[1][0].twinx()
#ax32.plot(x,TangentAlt, greycolour, linestyle = ':' )

ax[1][1].plot(dx[1:-1] , residual[1:-1] , 'b')
ax[1][1].set_title('Residual Doppler Shift (Hz)')
ax42 = ax[1][1].twinx()
ax42.plot(x,TangentAlt, greycolour, linestyle = ':' )
ax[1][1].grid(True)
ax42.set_label('Altitude of Tangent Point (km)')

#ax3 = fig.add_subplot(211)

plt.show()


print('stophere so the packages are still installed and u can test things')