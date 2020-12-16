import numpy as np
import spiceypy as spice
import spiceypy.utils.support_types as stypes
import pandas as pd
from os import *
import matplotlib.pyplot as plt
import time as timer
from scipy import constants
from tqdm import tqdm

import math

#Import custom modules
#import main
import atmosphere
import swiftmain

#THE CONTROL TIME = 636491202
#EXPERIEMNT mro->ODYSSEY  = 24269137.689745
#your egress value because you are an idiot = 241895765.6228938
def ElectricCall(time):
    process_id = getpid()
    print("Process ID:", process_id)
    epoch = 636491202.20059
    target = '-143'
    obs = '-41'
    #need a funtion that can interpolate between times at 10 Hz, then just sub those positions into the next function
    #MEX TGO positions to be created in this function with 10 Hz interpolation
    samplingfrequency =1
    referencedirection = [0,0,0]

    result = np.zeros([samplingfrequency])
    overshoot = np.zeros([samplingfrequency])
#636491202.20059
    #MEX,TGO, Distance = ephemerides(636491202,time, samplingfrequency)
    #mex = MEX[:,0] ; tgo = TGO[0]
    for i in range(samplingfrequency):
        [tgo1, _] = spice.spkpos('MARS', epoch - time, 'IAU_MARS', 'NONE', target) 
        [mex1, _] = spice.spkpos('MARS', epoch - time, 'IAU_MARS', 'NONE', obs)   
        dis = mex1-tgo1
        Distance = np.linalg.norm(dis)
        #print(f"Progess:{(i/samplingfrequency)*100} %")
        #mex = MEX[:,i]+0 ; tgo = TGO[:,i]+0
        mex = mex1 ; tgo = tgo1
        initialangle,xyzpoints= producegeometrymeter(mex,tgo) #make 2D for huge speed improvements.
        Bending, S , referencedirection = flatbending(xyzpoints,initialangle, mex,tgo, referencedirection)
        result = S
        #result[i] = np.stack((ElectricDistance, Distance), axis = 0) 
    return S

def TangentPointAltitude(time):
    epoch = 636491202.20059
            
    target = '-143'
    obs = '-41'
    alt = np.zeros(len(time)+1)
    for i in time:
        [tgo, _] = spice.spkpos('MARS', epoch - i, 'IAU_MARS', 'NONE', target) 
        [mex, _] = spice.spkpos('MARS', epoch - i, 'IAU_MARS', 'NONE', obs)  
        [states,_] = spice.spkezr(target, epoch-i, 'IAU_MARS', 'NONE', obs)
        sc2scvector = states[0:3]
        displacement = np.linalg.norm(sc2scvector)
        sc2scunitvector = np.true_divide(sc2scvector, displacement)
        marsrad = spice.bodvrd('MARS', 'RADII', 3)
        _,alt[i] = spice.npedln(marsrad[1][0], marsrad[1][1], marsrad[1][2], tgo, sc2scunitvector)
    return alt *1000

def GeoCall(time):
    epoch = 636491202.20059
    target = '-143'
    obs = '-41'
    process_id = getpid()
    print("Process ID:", process_id)
    [tgo1, _] = spice.spkpos('MARS', epoch - time, 'IAU_MARS', 'NONE', target) 
    [mex1, _] = spice.spkpos('MARS', epoch - time, 'IAU_MARS', 'NONE', obs)   
    dis = mex1-tgo1
    result = np.linalg.norm(dis)

    return result


#the smallest time quantum in SPICE is 1 second. to measure a dopplershift we need faster than this. 
# To sidestep this limmitation we interpolate 10 positions between seconds
def ephemerides(et,when, samplingfrequency):
    Distance = np.zeros([samplingfrequency])
    #Find the locations of MEX & TGO at epoch and epoch +1 second
    time = et+when # here is when u have set the order of the 
    TGO = np.zeros((3,samplingfrequency)) ; MEX = np.zeros((3,samplingfrequency))
    [tgo1, _] = spice.spkpos('MARS', time-samplingfrequency, 'IAU_MARS', 'NONE', target) ; [tgo2, _] = spice.spkpos('MARS', time, 'IAU_MARS', 'NONE', target)
    [mex1, _] = spice.spkpos('MARS', time-samplingfrequency, 'IAU_MARS', 'NONE', '-41')  ; [mex2, _] = spice.spkpos('MARS', time, 'IAU_MARS', 'NONE', '-41') 
    # [dis, _] = spice.spkpos('-143', time+1, 'IAU_MARS', 'NONE', '-41') 
    # Distance = np.linalg.norm(dis)
    #find ten positions between the two epochs for both MEX & TGO
    delta_tgo = (tgo2 - tgo1)/samplingfrequency ; delta_mex = (mex2-mex1) /samplingfrequency
    for i in range(samplingfrequency): 
        MEX[:,i] = mex1 + (delta_mex *i) 
        TGO[:,i] = tgo1 + (delta_tgo *i)
        #[dis, _] = spice.spkpos('-143', time+(i/samplingfrequency), 'IAU_MARS', 'NONE', '-41') 
        dis = MEX[:,i]-TGO[:,i]
        Distance[i] = np.linalg.norm(dis)
    return MEX, TGO, Distance





def producegeometrymeter(MEX,TGO):
#maybe completly thin this out, you know this is moslty pointless, what does it actually make

    class SpiceVariables:
        obs  = '-41' # NAIF code for MEX '-74'
        target = '-143'# NAIF code for TGO ['EARTH'/'SUN'/ a groundstation etc] 'MARS ODYSSEY'
        obsfrm = 'IAU_MARS'
        abcorr = 'NONE'
        crdsys = 'LATITUDINAL'
        coord  = 'LATITUDE'
        stepsz = 1.0 # Check every [300] seconds if there is an occultation
        MAXILV = 100000 #Max number of occultations that can be returned by gfoclt
        bshape = 'POINT'
        fshape = 'DSK/UNPRIORITIZED'
        front = 'MARS'
        fframe = 'IAU_MARS'
        TFMT = 'YYYY-MM-DD HR:MN:SC' # Format that Cosmographia understands

    sv = SpiceVariables()
    
    #THIS COULD BE REMOVED
            # [TGO, _] = spice.spkpos(sv.front, et-when, sv.fframe, 'NONE', sv.target)
            # [MEX, _] = spice.spkpos(sv.front, et-when, sv.fframe, 'NONE', sv.obs)
    TGO = TGO +0 ; MEX = MEX +0 #force to be non-strided
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
    #[SUN, _] = spice.spkpos(sv.front, et, sv.fframe, 'NONE', 'SUN')
    for i in range(dist):
        xyzpoint = MEX + (i * unitsc2sc) #move along ray, 1000 wavelength distance at a time (685 m). but unitsc2sc is in km...
        xyzpoints[:,i] = xyzpoint
        #sza[0,i] = spice.vsep(SUN,xyzpoint)
        angleprogression[0,i] = (spice.vsep( xyzpoint, MEX)) * (180 / math.pi)
        points[:,i] = spice.recgeo(xyzpoint, equatorialradii,flatteningcoefficient)
        points[0,i] = (points[0,i] * (-180 / math.pi))
        points[1,i] = (points[1,i] * (-180 / math.pi))
        
        

    # ray = np.concatenate((points,sza), axis=0) # important for when sza is included

    #plt.plot(angleprogression[0,:], ray[2,:])
    #plt.show()

    # ray is in lat/lon/alt + sza and xyzpoints is cartesian, both describe the same thing
    return  initialangle,xyzpoints




def flatbending(xyzpoints,initialangle, MEX,TGO, referencedirection):

    class SpiceVariables:
        obs  = '-74' # NAIF code for MEX
        target = 'MARS ODYSSEY'# NAIF code for TGO ['EARTH'/'SUN'/ a groundstation etc]
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

    sv = SpiceVariables()


    #form a coordinate system where tgo is @ y=0 and x= (5000 +norm), Mar's Barrycenter being @ [5000,0]
    subgroupsize = 1
    #initialise non-global variables
    miniray = np.zeros(subgroupsize)
    raystep = np.zeros((2,100000000))# create a large array to populate and then shrink later

    barry2mex = np.linalg.norm(MEX)
    barry2tgo = np.linalg.norm(TGO)

    #find the martian geomoerty so you can reliably find the altitude of a point
    marsrad = spice.bodvrd(sv.front, 'RADII', 3)
    flatteningcoefficient = ( marsrad[1][0] - marsrad[1][2] ) / marsrad[1][0]
    equatorialradii = marsrad[1][0]


    TGO = TGO +0 ; MEX = MEX +0 #force to be non-strided
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

    #THIS NEEDS TO VARY IF THERE IS AN OVERSHOOT
    unit = 1 # in km
    unitoriginal = unit
    
    rotationvector = np.array(( (np.cos(initialtheta), -np.sin(initialtheta)),
                                (np.sin(initialtheta),  np.cos(initialtheta)) ))

    #get unit vecotr of -MEX (then add this vecotr to MEX for each alt calcultation)
    unitmex = -mex/barry2mex #unit direction (2d)
    
    if referencedirection == [0,0,0]:
        initialdirection = unitmex.dot(rotationvector) * unit #make a 2d vector coming from MEX YOU DO NOT KNOW WHAT WAY THIS IS ROTATING
    else:# if there is a value for the fed-back starting direction than use this as the first firection
        initialdirection = referencedirection

    iterationcount =0
    #while iterationcount<100:
    #print( "Finding Bending Angle (", str(iterationcount) ,"% Complete)")
    errorstore = np.zeros((11,100000))
    S = np.zeros(20000)
    
    #IF REFERCEDRECTION==0 DO NORMAL, IF /= INCLUDE THIS AS THE FIRST DIRECTION.

    while iterationcount <10:
        tic  = timer.perf_counter()
        
        if iterationcount == 0:
            direction = initialdirection
            
        else:
            
            missangle = missangle/1

            missrotationvector = np.array(( (np.cos(missangle), -np.sin(missangle)),
                                            (np.sin(missangle),  np.cos(missangle)) ))

            direction = initialdirection.dot(missrotationvector)
            #check the differecne between the initial and the direction, see if the same for both
            CHECK_ME = direction - initialdirection
            initialdirection = direction 

        turningcounter=0
        stage =0
        t=0
        unit = unitoriginal
        #with tqdm(total = mex[1], desc = "Progress", leave=False) as pbar:
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



                #IF THE Y VALUE DROPS BELOW 0, LOOP BACK WITH A STEP SIZE OF 1/10TH#################################################################<- HERE
                if point[1] < 0: #if the position drops below the x axis
                    direction = direction/10 #drop the step size
                    unit= unit/10
                    #stage = stage+1

                    #MAYBE SHRINK THE DIRECTION INSTEAD OF THE UNIT, IT DOESNT GET REINITIALED IN THIS WHILE LOOP
                    # t is not incrememented so that it can loop back to the last position 
                    # , t-1 is the position just before crossing over into negative space 
                    if abs(point[1]) < 0.00001: # is it smaller than cm
                        stage = stage+1 #increase the stage value so the while loop is exited
                    continue
                        
                    



                        # #this section allows for better timing of the function, increment the progresbar by 1 if the current 
                        # #position goes one y-value lower
                        # currenty = raystep[1,t]
                        # progress[1] = mex[1] - currenty #this value will be increasing from 0 -> Mex height
                        # increment = np.floor(progress[1])- np.floor(progress[0]) #only when 
                        # if increment ==1 :
                        #     pbar.update(1)
                        # progress[0] = progress[1]


                if abs(N1) < 1e-20:#catch for precision errors
                    S[t] = unit - ( N1 * unit) # whilst neglegible N, Electric distance is can simply be inversly proportional to N
                    t=t+1
                    N0=N1
                    continue
                
                #THE SECTION BELOW IS ONLY ACCESSED WHEN THE REFRACTIVTY IS ABOVE E-20
                #print('Current Y possition is', currenty) #this is alt, so is ~3389 km smaller than the vector
                r= N0/N1 #NEED MORE PRECISION
                numorator = N0+1
                denominator = N1+1
                rbending = numorator/denominator #average would just add 1 to total N
                # if t==5000: #only bend when there is a refractive gradient between consecutive air volumes[NO CHANGE IN DIRECTION]
                #     t=t+1
                #     N0=N1
                #     continue

                #this section is only reached if a turning is going to happen
                # ,testing to see if there is 10 X less turing if the units are 10X smaller -TRUE
                TEST = float(rbending)
                if TEST != 1:
                    turningcounter = turningcounter+1


                # !! NOW WITH PRECISION !!

                #find the angle between the unit (air volume) boarder and the current direction
                unitrotationaxis = raystep[:,t]/((np.linalg.norm(raystep[:,t]))*unit)
                #unitdirection = direction #MAYBE ALTERING UNITS WILL EFFECT THIS * OR / BY UNIT, CANT FIGURE OUT NOW, OK WHEN UNIT =1
                DotProduct=  np.dot((unitrotationaxis),direction)
                AngleofIncidence = (math.pi/2)- np.arccos(DotProduct) #angle it enters the next air volume
                #simple snell law to find the bending angle (should be tiny angle)     
                AngleofRefraction = np.arcsin(rbending * np.sin(AngleofIncidence))
                # THIS IS NOT EXACTLY WHAT THE TURN IN DIRECTION IS, NEED TO THINK ABOUT
                rotateby = ((AngleofIncidence - AngleofRefraction))#+ve =clockwise, -ve=anticlockwise

                INCIDENCEDEGREES = np.degrees(AngleofIncidence)
                REFRACTIONDEGREES = np.degrees(AngleofRefraction)
                ROTATIONDEGREES = np.degrees(rotateby)

                #an if statement is required, if this 
                if ROTATIONDEGREES>1 or ROTATIONDEGREES<-1:
                    print('stophere, u r bending to much')
                rotationvector = np.array(( (np.cos(rotateby), -np.sin(rotateby)),
                                            (np.sin(rotateby),  np.cos(rotateby))))
                direction = direction.dot(rotationvector)

                
                N0=N1

                #store N1 to calc electric distance 
                S[t] = unit - ( N1 * unit)
                t=t+1

            #pbar.refresh()
        
        unit_initial = initialdirection/ np.linalg.norm(initialdirection)
        dot_product = np.dot(unit_initial,unitmex)
        FinalBendingAngle = np.arccos(dot_product)

        error = np.zeros(t)   
        #print("Number of turns:", turningcounter)
        #error = swiftmain.finderror(raystep, UnrefractedRay) # also find the y-overshoot here

        miss = error # 1D along the abscissa

        #update for the miss angle, going to include miss in the ordinate
        #START EDITING HERE
        miss = point[0] - UnrefractedRay[0,-1]# X domain

        deltaX = UnrefractedRay[0,-1] #TGO X value
        deltaY = UnrefractedRay[1,0] # this is the height of the whole 2d scene (both refract and unrefacted have the same height)
        unrefractedangle = np.arctan(deltaX/deltaY)
        refractedangle = np.arctan((deltaX + miss)/deltaY)
        missangle = refractedangle - unrefractedangle # if positive, then rotate clockwise
        #missangle = (np.arcsin(miss/UnrefractedDistance)) # this shouldnt work
        toc  = timer.perf_counter()
        passingtime = toc-tic
        #print('miss =', format(miss*1000, '.5f') ,'m || Angle =', np.degrees(missangle) ,
        #                    '° || Speed =',passingtime,' Sec \n', sep = " ", end= " ", flush =True)
        
        if abs(miss) < 1e-4: #is the miss smaller than 10 cm?
            #ploterrortraces(errorstore,t)
            break

        iterationcount = iterationcount +1 

    #find the total bending angle at MEX for this final configuration


    #from the refractive profile from MEX ->TGO, calc the intergral of the change in wavelength to aquire doppler
    S =  S[S != 0]
    ElectricDistance = (np.sum(S))#N * wavelengths in a km (UNITS MUST BE KEPT THE SAME AS STEP-SIZE)[+ OVERSHOT BECAUSE IT IS A NEGATIVE VARIABLE ]

    return FinalBendingAngle, ElectricDistance , initialdirection #feedback the starting vector for speed



#NOW ONLY PRODUCING THE ERROR ON THE FINAL VALUE
def finderror(refract, unrefract):
    #the deviation from the straight line can be found by rotating the refracted propergation by the original starting
    # direction. then subtract the MEX x-value. this will produce the refracted ray deviating from a verticle straight line. 
    # This method removes any rounding or index searching errors
    
    #find angle of unrefracted (previous anlgle not accurate enough found in 3d so slightly inaccurate)
    unrefracted_start2end = unrefract[:,0] - unrefract[:,-1] 
    unrefract_end  = unrefract[:,-1]
    unrefracted_start = [0,unrefract[1,0]] #this should now be a verticle line
    unit_start2end = unrefracted_start2end/np.linalg.norm(unrefracted_start2end)
    unit_start = unrefracted_start/np.linalg.norm(unrefracted_start)
    dot_product = np.dot(unit_start2end,unit_start)
    newangle = np.arccos(dot_product)# should be smaller

    #refract must be shortened to its no 0 components (BAD METHOD DOING THIS)
    refractx= refract[0,:] ; refracty= refract[1,:]
    refractx = refractx[refractx !=0] ; refracty = refracty[refracty !=0]
    refract = np.vstack((refractx,refracty))

    yovershoot = refracty[-1] #should be -ve
    uwot = unrefract[0,-1]#NO BENDING
    uwot2 = refractx[-1]
    xovershoot = unrefract[0,-1] - refractx[-1]# should be +ve [net bending upwards with dominant iono]
    #THIS NEEDS TO BE A HYPO
    overshoot = np.sqrt((refracty[-1]**2) + ((unrefract[0,-1] - refractx[-1])**2))# sqrt(y_overshoot^2 + x_overshoot^2)

    #or the simple version to see if it is the hypno function that is screwed
    #overshoot = refracty[-1]

    ox, oy = unrefract[:,0] #extract MEX as origin
    px, py = refract[:,-1]#in swiftmain, we only want the final value
    rotx = ox + math.cos(-newangle) * (px - ox) - math.sin(-newangle) * (py - oy)#rotate each point by initial starting angle, so it goes downwards
    error = rotx

    error = error - unrefract[0,0] #subtract the origin x-value from every value in error

    return error

# Find the average refractivty of volume intercepted by ray
def findrefractivity(alt, step): 
    #must convert these to 32decimal
    ionoresidual = atmosphere.iono(alt,step)
    #ionoresidual = 0
    neutralresidual = 0
    #neutralresidual = atmosphere.neutral(alt,step)
    n = np.sum( (ionoresidual + neutralresidual))/step #produce the integral refractive index and divide by lenght to find average
    return n

# this section is required to find the occulation epoch
def core():

    class SpiceVariables:
        obs  = '-74' # NAIF code for MEX
        target = 'MARS ODYSSEY'# NAIF code for TGO ['EARTH'/'SUN'/ a groundstation etc]
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

    sv = SpiceVariables()

    #-----------------------------------------------------<VALUES TO EDIT REGULARLY>----------------------------------------
    # If you only wish to analysis mutual [cross-link] occultation between MEX and TGO, then this is the only section that
    # needs to be edited
    start  = '2020 MAR 1'
    stop   = '2020 MAR 3'
    OCCSELECTION = 17 # Which occultation do you wish to see in Cosmographia? [optional]
    here = path.abspath(path.dirname(__file__))
    PathtoMetaKernel1 = here + '/TGO/krns/mk/em16_plan.tm'
    PathtoMetaKernel2 = here + '/MEX/krns/mk/MEX_OPS.tm'
    #-----------------------------------------------------------------------------------------------------------------------

    spice.furnsh(PathtoMetaKernel1)
    spice.furnsh(PathtoMetaKernel2)

    sv = SpiceVariables()

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



    occs = pd.DataFrame(ingresslist, columns=['Time'])
    occ = occs.Time[OCCSELECTION]

    return occ

# sv = main.SpiceVariables()
# occ = core()
#print("strange result:", occ)














                                        # #print(result)





                                        # #form the dataframe
                                        # length = 120
                                        # occs = pd.DataFrame(ingresslist, columns=['Time'])
                                        # residualdoppler = np.zeros(length)
                                        # velocitydoppler = np.zeros(length)
                                        # RESIDUALSUM = np.zeros(length)
                                        # #Calculate the residual doppler as the sum of the neutral and ionosphere
                                        # tic  = timer.perf_counter()
                                        # for time in tqdm(range(length)): #begin time at occultation epoch and go to 2 mins pior


                                        #     ray, dist, totalperiods, remainingdistance = main.producegeometrylamda(occs.Time[OCCSELECTION], sv, time*8)# Produce geometry in wavelenghts to investigate electric distance

                                        #     _,_,_,_, alt = main.Location(occs.Time[OCCSELECTION], sv, time*8)

                                        #     ionoresidual = atmosphere.iono(ray[2,:],totalperiods)
                                        #     neutralresidual = atmosphere.neutral(ray[2,:],totalperiods)
                                        #     residual = 1 + (ionoresidual + neutralresidual)

                                        #     # plt.plot(range(totalperiods),residual[0,:])
                                        #     # plt.title("Refractive Index through Propergation of Ray")
                                        #     # plt.xlabel("MEX->TGO distance (km)")
                                        #     # plt.ylabel("Refractive Index")
                                        #     # plt.show()

                                        #     #DO A TEST TO SEE IF THE NET REFRACTIVE INDEX CHANGES OVER TIME, THEN U CAN HONE IN ON THE ERRRO
                                        #     # account for the plus 1 

                                        #     [electricdistance, geometric, resdopplershift] = main.doppler(residual, totalperiods, dist, remainingdistance)
                                        #     miss = electricdistance - dist + remainingdistance # a possitive number as electric distance is greater that geometric due to iono
                                        #     howwrongurcodeis = geometric - dist # is this purly due to rounding of that 9945 (each 1 is ~685 m)
                                        #     residualdoppler[time] = resdopplershift

                                            

                                        #     #find the geometric doppler too UNTESTED MODULE
                                        #     sc2scstates = spice.spkezr(sv.target, (occs.Time[OCCSELECTION] - time*8), sv.fframe, 'LT+S', sv.obs)
                                        #     velocityvector = sc2scstates[0][3:6]
                                        #     velocityvector = velocityvector[0:3]
                                        #     positionalvector =  sc2scstates[0][0:3]
                                        #     positionalvector = positionalvector[0:3]
                                        #     velocityangle = spice.vsep( positionalvector, velocityvector) #rads
                                        #     relativevelocity = np.linalg.norm(velocityvector) * np.cos(velocityangle) 
                                        #     geometricdopplershift = -(relativevelocity/constants.c) * 437.1e6
                                        #     velocitydoppler[time] = geometricdopplershift *1000 # becuase spice is in km and c is in m



                                        # toc = timer.perf_counter()
                                        # passingtime = toc-tic
                                        # print('elapsed time in seconds:', passingtime)
                                        # noise = np.random.normal(0,50,velocitydoppler.shape)
                                        # RESIDUALSUM = residualdoppler + velocitydoppler + noise

                                        # fig, ax = plt.subplots(3)
                                        # ax[0].plot(range(-960,0,8),velocitydoppler[: : -1] )
                                        # ax[0].set_title('Geometric', loc='left')
                                        # ax[0].set_xlabel('Time after Occ (s)')
                                        # ax[0].set_ylabel('Doppler Shift (Hz)')
                                        # ax[1].plot(range(-960,0,8),residualdoppler[: : -1] )
                                        # ax[1].set_title('Residual', loc='left')
                                        # ax[1].set_xlabel('Time after Occ (s)')
                                        # ax[1].set_ylabel('Doppler Shift (Hz)')
                                        # ax[2].plot(range(-960,0,8),RESIDUALSUM[: : -1])
                                        # ax[2].set_title('Product + Noise', loc='left')
                                        # ax[2].set_xlabel('Time to Horizon Epoch (s)')
                                        # ax[2].set_ylabel('Doppler Shift (Hz)')
                                        # plt.show()

                                        # #then add noise