import numpy as np
import spiceypy as spice
from scipy import constants
from scipy.interpolate import InterpolatedUnivariateSpline
import math
import time as timer

#from decimal import *
from tqdm import tqdm
from mpmath import *
mp.dps = 10 ; mp.pretty = True 
#getcontext().prec = 35 #set the precision for the bending in function "flatbending" 



import JSONtemplate as JS
import pandas as pd
import json
from PIL import Image
import matplotlib.pyplot as plt
import cartopy.crs as ccrs #import the coordinate refernece system
from scipy.spatial.transform import Rotation as R

import atmosphere


class SpiceVariables:
    obs  = '-41' # NAIF code for MEX
    target = '-143'# NAIF code for TGO ['EARTH'/'SUN'/ a groundstation etc]
    obsfrm = 'IAU_MARS'
    abcorr = 'NONE'
    crdsys = 'LATITUDINAL'
    coord  = 'LATITUDE'
    stepsz = 100.0 # Check every 300 seconds if there is an occultation
    MAXILV = 100000 #Max number of occultations that can be returned by gfoclt
    bshape = 'POINT'
    fshape = 'DSK/UNPRIORITIZED'
    front = 'MARS'
    fframe = 'IAU_MARS'
    TFMT = 'YYYY-MM-DD HR:MN:SC' # Format that Cosmographia understands



def producegeometrymeter(et,sv,when):
    [TGO, _] = spice.spkpos(sv.front, et-when, sv.fframe, 'NONE', sv.target)
    [MEX, _] = spice.spkpos(sv.front, et-when, sv.fframe, 'NONE', sv.obs)

    dist = math.floor(spice.vdist(TGO,MEX))

    angleseparation = (spice.vsep(MEX, TGO)) # angle taken a mars center
    initialangle  = (spice.vsep(-MEX, (TGO-MEX))) * (180/math.pi)# angle taken at mars-MEX-tgo, that points to tgo. needed for the bending functions original starting angle  
    #script needs to work via periods of ray and not meters. [totalperiods is the main iterable, not meters]
    
    sc2sc = TGO - MEX
    norm = np.linalg.norm(sc2sc)
    unitsc2sc = sc2sc/norm #this needs to shrink if the repeatable expands
    points = np.empty([3,dist])
    sza = np.empty([1,dist])
    angleprogression = np.empty([1,dist])

    xyzpoints = np.zeros([3,dist])
    marsrad = spice.bodvrd(sv.front, 'RADII', 3)
    flatteningcoefficient = ( marsrad[1][0] - marsrad[1][2] ) / marsrad[1][0]
    equatorialradii = marsrad[1][0]
    # find direction of sun, it will not change much during the occultation. so only calc it once
    [SUN, _] = spice.spkpos(sv.front, et, sv.fframe, 'NONE', 'SUN')
    for i in range(dist):
        xyzpoint = MEX + (i * unitsc2sc) #move along ray, 1000 wavelength distance at a time (685 m). but unitsc2sc is in km...
        xyzpoints[:,i] = xyzpoint
        sza[0,i] = spice.vsep(SUN,xyzpoint)
        angleprogression[0,i] = (spice.vsep( xyzpoint, MEX)) * (180 / math.pi)
        points[:,i] = spice.recgeo(xyzpoint, equatorialradii,flatteningcoefficient)
        points[0,i] = (points[0,i] * (-180 / math.pi))
        points[1,i] = (points[1,i] * (-180 / math.pi))
       
        
        

    ray = np.concatenate((points,sza), axis=0) # important for when sza is included

    #plt.plot(angleprogression[0,:], ray[2,:])
    #plt.show()

    # ray is in lat/lon/alt + sza and xyzpoints is cartesian, both describe the same thing
    return  initialangle, MEX,TGO,  xyzpoints



def flatbending(xyzpoints,initialangle, sv, MEX,TGO):
    #form a coordinate system where tgo is @ y=0 and x= (5000 +norm), Mar's Barrycenter being @ [5000,0]
    subgroupsize = 10
    #initialise non-global variables
    miniray = np.zeros(subgroupsize)
    raystep = np.zeros((2,100000000))# create a large array to populate and then shrink later

    barry2mex = np.linalg.norm(MEX)
    barry2tgo = np.linalg.norm(TGO)

    #find the martian geomoerty so you can reliably find the altitude of a point
    marsrad = spice.bodvrd(sv.front, 'RADII', 3)
    flatteningcoefficient = ( marsrad[1][0] - marsrad[1][2] ) / marsrad[1][0]
    equatorialradii = marsrad[1][0]

    _,_, MEXalt = spice.recgeo(MEX, equatorialradii,flatteningcoefficient)
    _,_, TGOalt = spice.recgeo(TGO, equatorialradii,flatteningcoefficient)

    #the possition of MEX is found by assuming that it will be somewhere over the relative horizon from TGO 
    # (meaning over θ = 90°), finding the angle between the MEX and TGO's negative vector, will give the coords of MEX
    MexRelativeElevation = spice.vsep(-TGO, MEX) #radians
    mex_y = barry2mex * np.sin(MexRelativeElevation)
    mex_x = barry2mex * np.cos(MexRelativeElevation)

    mex = np.array([0-mex_x, mex_y])
    tgo = np.array([0+barry2tgo, 0])
    barry = np.array([0, 0])

    #to plot the non-refracted propogation, we must convert the 3d xyzpoints to 2d, we do this the same way we found the x&y for MEX 
    # ,using the norm distance and sep from -TGO5
    length = np.size(xyzpoints,1)
    UnrefractedDistance = np.linalg.norm(xyzpoints[:,0]- xyzpoints[:,-1]) #in km
    UnrefractedRay = np.zeros([2,length])
    for i in range(length): #conversion to 2D
        point = xyzpoints[:,i] + 0 #need to put vector into temp variable as spice cant handle strided array inputs
        angle = spice.vsep(-TGO, point)
        norm = np.linalg.norm(point)
        point_x = norm * np.cos(angle)
        point_y = norm * np.sin(angle)
        UnrefractedRay[0,i] = 0 - point_x
        UnrefractedRay[1,i] = point_y

     #this will produce and angle that is likly not going to be exactly on 
    #the original propagation path, you compare to this if there is a drifting error, as both this and the resultant refracted ray 
    # have the same bias error. THIS ANGLE IS BENDING ANTICLOCKWISE IN THIS FRAME (BENDING UPWARDS)
    initialtheta = -(spice.vsep(MEX-TGO, MEX))
    nicetohave = np.degrees(initialtheta)

    unit = 1 # in km
    
    
    rotationvector = np.array(( (np.cos(initialtheta), -np.sin(initialtheta)),
                                (np.sin(initialtheta),  np.cos(initialtheta)) ))

    #get unit vecotr of -MEX (then add this vecotr to MEX for each alt calcultation)
    unitmex = -mex/barry2mex #unit direction (2d)
    initialdirection = unitmex.dot(rotationvector) * unit #make a 2d vector coming from MEX YOU DO NOT KNOW WHAT WAY THIS IS ROTATING


    iterationcount =0
    #while iterationcount<100:
    #print( "Finding Bending Angle (", str(iterationcount) ,"% Complete)")
    errorstore = np.zeros((11,100000))
    Nstore = np.zeros(20000)
    
    while iterationcount <100:
        tic  = timer.perf_counter()
        if iterationcount ==0:
            direction = initialdirection
        else:
            # the initail direction must be rotated by a 10th of the miss at the end
            missangle = missangle/1
            missrotationvector = np.array(( (mp.cos(missangle), mp.sin(missangle)),
                                            (mp.sin(missangle),  mp.cos(missangle)) ))
            direction = initialdirection.dot(missrotationvector)
            initialdirection = direction 

        turningcounter=0
        progress= [0,0]
        stage =0
        t=0
        with tqdm(total = mex[1], desc = "Progress", leave=False) as pbar:
            while stage < 2: #==0first unit, so move two units. ==1 propergate step by step. ==2 exit and analyse entire path 
                #take the alt at 10 possitions across this unit

                #lets get a quick calculated for the magnitude of the direction
                MAAAAG = np.linalg.norm(direction) # this should = unit

                if stage==0:
                    for k in range(subgroupsize): #this starts with 0
                        point = mex + ((k+1) *(direction/subgroupsize))
                        #_,_,miniray[k] = spice.recgeo(point, equatorialradii,flatteningcoefficient)
                        miniray[k] =  np.linalg.norm(point) - 3389 #average radii of mars
                    N0 = findrefractivity(miniray,subgroupsize)
                    raystep[:,t] = point #save the last location
                    t=t+1
                    stage =stage+1

                if stage==1:
                    for k in range(subgroupsize):
                        point = raystep[:,t-1] + ((k+1) *(direction/subgroupsize)) #am i double counting the end of the last and the start of the next?
                        #_,_,miniray[k] = spice.recgeo(point, equatorialradii,flatteningcoefficient)  #THIS ONLY WORKS IN 3D
                        # IMPLEMENTING MARS AS A SIMPLE CIRCLE OF AVERAGE 3389 KM RADIUS, !THIS WILL BE UPDATED TO ELLIPSE!
                        miniray[k] =  np.linalg.norm(point) - 3389
                    raystep[:,t] = point#9 is the end of the unit, and refraction always happens relative to the center of refractivity, so rotate off this vector
                    N1 = findrefractivity(miniray,subgroupsize)

                    if point[1] < 0: #if the position drops below the x axis
                        stage = stage+1 #increase the stage value so the while loop is exited

                    #this section allows for better timing of the function, increment the progresbar by 1 if the current 
                    #position goes one y-value lower
                    currenty = raystep[1,t]
                    progress[1] = mex[1] - currenty #this value will be increasing from 0 -> Mex height
                    increment = np.floor(progress[1])- np.floor(progress[0]) #only when 
                    if increment ==1 :
                        pbar.update(1)
                    progress[0] = progress[1]


                    if abs(N0) < 1e-20:#catch for precision errors
                        Nstore[t] = N1
                        t=t+1
                        N0=N1
                        continue
                    
                    #print('Current Y possition is', currenty) #this is alt, so is ~3389 km smaller than the vector
                    r= N0/N1 #NEED MORE PRECISION
                    numorator = mpf(N0)+mpf(1)
                    denominator = mpf(N1)+mpf(1)
                    rbending = mp.fdiv(numorator,denominator) #average would just add 1 to total N
                    if t==5000: #only bend when there is a refractive gradient between consecutive air volumes[NO CHANGE IN DIRECTION]
                        t=t+1
                        N0=N1
                        continue

                    #this section is only reached if a turning is going to happen
                    # ,testing to see if there is 10 X less turing if the units are 10X smaller -TRUE
                    TEST = float(rbending)
                    if TEST != 1:
                        turningcounter = turningcounter+1


                    # !! NOW WITH PRECISION !!

                    #find the angle between the unit (air volume) boarder and the current direction
                    unitrotationaxis = raystep[:,t]/np.linalg.norm(raystep[:,t])
                    #unitdirection = direction #MAYBE ALTERING UNITS WILL EFFECT THIS * OR / BY UNIT, CANT FIGURE OUT NOW, OK WHEN UNIT =1
                    DotProduct=  np.dot(unitrotationaxis,direction)
                    AngleofIncidence = (math.pi/2)- mp.acos(DotProduct) #angle it enters the next air volume
                    #simple snell law to find the bending angle (should be tiny angle)     
                    AngleofRefraction = mp.asin(rbending * mp.sin(AngleofIncidence))
                    # THIS IS NOT EXACTLY WHAT THE TURN IN DIRECTION IS, NEED TO THINK ABOUT
                    rotateby = ((AngleofIncidence - AngleofRefraction))#+ve =clockwise, -ve=anticlockwise

                    INCIDENCEDEGREES = mp.degrees(AngleofIncidence)
                    REFRACTIONDEGREES = mp.degrees(AngleofRefraction)
                    ROTATIONDEGREES = mp.degrees(rotateby)

                    #an if statement is required, if this 
                    if ROTATIONDEGREES>1 or ROTATIONDEGREES<-1:
                        print('stophere, u r bending to much')
                    rotationvector = np.array(( (mp.cos(rotateby), -mp.sin(rotateby)),
                                                (mp.sin(rotateby),  mp.cos(rotateby))))
                    direction = direction.dot(rotationvector)

                    
                    N0=N1

                    #store N1 to calc electric distance 
                    Nstore[t] = N1
                    t=t+1

            #pbar.refresh()

        error = np.zeros(t)   
        #print("Number of turns:", turningcounter)
        error = finderror(raystep, UnrefractedRay)
        miss = error

        missangle = mp.asin(miss/ UnrefractedDistance)
        toc  = timer.perf_counter()
        passingtime = toc-tic
        print(' miss =', format(miss*1000, '.5f') ,'m || Angle =', nstr(mp.degrees(missangle),5) ,
                            '° || Speed =',passingtime,' Sec \n', sep = " ", end= " ", flush =True)
        if abs(miss) <1e-5: #is the miss smaller than 10 cm?
            #ploterrortraces(errorstore,t)
            break
        iterationcount=iterationcount+1

    #find the total bending angle at MEX for this final configuration
    unit_initial = initialdirection/ np.linalg.norm(initialdirection)
    dot_product = np.dot(unit_initial,unitmex)
    FinalBendingAngle = mp.acos(dot_product)

    #from the refractive profile from MEX ->TGO, calc the intergral of the change in wavelength to aquire doppler
    Nstore =  Nstore[Nstore != 0]
    Nstore = Nstore + 1 #convert N deltas into values of N
    ElectricDistance = np.sum(Nstore) #N * wavelengths in a km (UNITS MUST BE KEPT THE SAME AS STEP-SIZE)

    return FinalBendingAngle, ElectricDistance



#NOW ONLY PRODUCING THE ERROR ON THE FINAL VALUE
def finderror(refract, unrefract):
    #the deviation from the straight line can be found by rotating the refracted propergation by the original starting
    # direction. then subtract the MEX x-value. this will produce the refracted ray deviating from a verticle straight line. 
    # This method removes any rounding or index searching errors
    
    #find angle of unrefracted (previous anlgle not accurate enough found in 3d so slightly inaccurate)
    unrefracted_start2end = unrefract[:,0] - unrefract[:,-1] 
    unrefracted_start = [0,unrefract[1,0]] #this should now be a verticle line
    unit_start2end = unrefracted_start2end/np.linalg.norm(unrefracted_start2end)
    unit_start = unrefracted_start/np.linalg.norm(unrefracted_start)
    dot_product = np.dot(unit_start2end,unit_start)
    newangle = np.arccos(dot_product)# should be smaller

    #refract must be shortened to its no 0 components (BAD METHOD DOING THIS)
    refractx= refract[0,:] ; refracty= refract[1,:]
    refractx = refractx[refractx !=0] ; refracty = refracty[refracty !=0]
    refract = np.vstack((refractx,refracty))

    ox, oy = unrefract[:,0] #extract MEX as origin
    px, py = refract[:,-1]#in swiftmain, we only want the final value
    rotx = ox + math.cos(-newangle) * (px - ox) - math.sin(-newangle) * (py - oy)#rotate each point by initial starting angle, so it goes downwards
    error = rotx

    error = error - unrefract[0,0] #subtract the origin x-value from every value in error

    return error

# Find the average refractivty of volume descriped by ray
def findrefractivity(alt, step): 
    #must convert these to 32decimal
    ionoresidual = atmosphere.iono(alt,step)
    neutralresidual = atmosphere.neutral(alt,step)
    n = np.sum( (ionoresidual + neutralresidual))/step #produce the integral refractive index and divide by lenght to find average
    return n