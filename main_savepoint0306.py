import numpy as np
import spiceypy as spice
from scipy import constants
from scipy.interpolate import interpolate
import math
import array as arr
import shapely
from decimal import *
from tqdm import tqdm_gui 

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



# OCCGRADIENT forms the shape of the occultation profile by creating and 3 column array of the profile's
# cartesian coordinated in the Martian reference frame[x y z].This functions works backwards from the
# point of occulatation. the time intervals are seconds. The top of the profile is defined
# when the TGO [target] becomes closer than the tangent point.
def occgradient(front, time, fframe, obs, target):

    marsrad = spice.bodvrd(front, 'RADII', 3)# Extract the triaxial dimensions of Mars
    increaser = 0
    altlist = 0

    for i in range(2000):
        # Find relative positions of TGO and MEX
        [targetpos, _] = spice.spkpos(front, (time - i), fframe, 'LT+S', target)
        [sc2scvector, _] = spice.spkpos(target, (time - i), fframe, 'NONE', obs)
        [obspos, _] = spice.spkpos(front, (time - i), fframe, 'LT+S', obs)

        # Find the unit vector between the SCs
        displacement = math.sqrt(((sc2scvector[0]) ** 2) + ((sc2scvector[1]) ** 2) + ((sc2scvector[2]) ** 2))
        unitvector = np.true_divide(sc2scvector, displacement)

        # Find the point this unit vector is closest to the Mars
        [profilesurfacepoint, alt] = spice.npedln(marsrad[1][0], marsrad[1][1], marsrad[1][2], targetpos, unitvector)
        altlist = np.append(altlist, alt)

        # Find distance between the tangent point and MEX (obs)
        tangent2mexdist = np.linalg.norm(obspos - profilesurfacepoint)

        # Need to add the altitude to the surface bound point. 'profilesurfacepoint' is a vector
        # from Mars' barrycenter to the surface, therefore
        # The vector is also the direction of altitude, so product must be added to the surface point
        tangentpointunitvector = profilesurfacepoint / np.linalg.norm(profilesurfacepoint)
        tangentpoint = (tangentpointunitvector * alt) + profilesurfacepoint

        #the SPICE function 'npedln' does not calculate along a vector between the two SCs, the vector continues through TGO. This is correced
        #by the following. If the tangent2mex distance > than the tgo2mex distance, then the vector has passed though the TGO. End function when this occurs
        if tangent2mexdist > (displacement-50):
            tangentpoint = targetpos
            highesttime = i
            break

        #The peak height can also be defined by when the profile begins to planteau [when the profile raises by less than
        # 50 m/s ]
        if altlist[i-1]+0.5 >altlist[i]:
            highesttime = i

        # Form array of profile coordinates
        if i == 0:
            x = int(tangentpoint[0])
            y = int(tangentpoint[1])
            z = int(tangentpoint[2])
        else:
            x = np.append(x, int(tangentpoint[0]))
            y = np.append(y, int(tangentpoint[1]))
            z = np.append(z, int(tangentpoint[2]))

        # If there is a particular altitude that needs to be highlighted, can be done like this (eg 200km)
        if alt > 200 + increaser:
            lowesttime = i
            increaser = 100000

    # Form a array of cartesian coordinates, transpose for convinience
    profile = [x.T, y.T, z.T]

    return profile, lowesttime, highesttime




def CosmographiaCatalogFormer(et,sv):

    # Select and occultation of interest, calculate the shape of it's profile, then add this to a dateframe

    result = occgradient(sv.front, et, sv.fframe, sv.obs, sv.target)
    profile = result[0]
    lowesttime = result[1]
    highesttime = result[2]
    lowesttimeFMT = spice.timout((et - result[1]), sv.TFMT)
    highesttimeFMT = spice.timout((et - result[2]), sv.TFMT)
    endtimeFMT = spice.timout(et, sv.TFMT)
    profile = np.transpose(profile)
    profiledataframe = pd.DataFrame(profile[:][0:highesttime], columns=['X', 'Y', 'Z'])
    # Example polygon velocities [this effects the dynamic texture of the profile, purely aesthetic]. These values have
    # been chosen to be slow and variable
    vs = np.array([[1.4236, -2.4657, 6.3948],
                   [1.6404, -2.2997, 6.4047],
                   [1.8386, -2.1150, 6.4145],
                   [2.0166, -1.9136, 6.4243],
                   [2.1730, -1.6977, 6.4339],
                   [2.3068, -1.4696, 6.4435],
                   [2.4170, -1.2315, 6.4530],
                   [2.5029, -0.9859, 6.4625],
                   [2.5029, -0.9859, 6.4625],
                   [2.5029, -0.9859, 6.4625]])

    # Iterate the polygon velocity states to be the length of profile and then combine with the positions
    vs = vs.repeat(200, axis=0)
    vt = vs[:][0:highesttime]
    velocitiesdataframe = pd.DataFrame(vt, columns=['dX', 'dY', 'dZ'])
    finalprofile = pd.concat([profiledataframe, velocitiesdataframe], axis=1)

    # Construct a JSON template depending on the size of the profile, split into sections of 10 points [smaller sections,
    # more smoothly the profile forms over time]
    blockcount = math.floor(highesttime / 10) * 10
    i = 0
    while i < (blockcount / 10):
        JS.JSONiterated = JS.JSONiterated + JS.JSONiterable
        i = i + 1
    JSONtemplate = JS.JSONstart + JS.JSONiterated + JS.JSONend

    template = json.loads(JSONtemplate)  # load the template created above  so it can be editted as JSON

    itemcounter = 0

    # convert states into single string for cosmographia to understand [cosmogrpahia computes as
    # 3 postions and 3 velocities, no need for '\n']
    for i in range(0, blockcount, 10):
        block = (finalprofile[i:i + 10].to_numpy()) * (-1)  # inverse due to reference frame inversion
        ProfileBlock = block.tolist()  # convert to list for the JSON format

        states = []
        for p in range(10):
            states.extend(
                ProfileBlock[p])  # turn states into one long line to remove the '\n', which confuses cosmographia

        # vary the spread and intensity of the profile as it increases with alitude and the atmosphere becomes less dense
        profilespread = math.floor(itemcounter * (5 / (blockcount / 10)))
        profileintensity = math.floor(itemcounter * (30 / (blockcount / 10)))

        # edit the large JSON file iterativly, with 'itemcounter'moving down the JSON items (blocks of 10 profile points).
        # adding the states(position and velocity of the polygons), name of item (so comsmographia doesnt overwrite them),
        # the start time, end time(epoch of occultation), the spread and finally the intensity of colour
        template['items'][itemcounter]['geometry']['emitters'][0]['generator']['states'] = states
        template['items'][itemcounter]['name'] = 'Profile ' + str(itemcounter)
        template['items'][itemcounter]['geometry']['emitters'][0]['startTime'] = spice.timout((et - i), sv.TFMT)
        template['items'][itemcounter]['geometry']['emitters'][0]['endTime'] = endtimeFMT
        template['items'][itemcounter]['geometry']['emitters'][0]['velocityVariation'] = (profilespread + 1)
        template['items'][itemcounter]['geometry']['emitters'][0]['endSize'] = profileintensity
        itemcounter = itemcounter + 1

    # serialise the formed profile into a .json that cosmographia can read
    with open(' Profile.json', 'w') as file:
        json.dump(template, file, indent=3)


# Find the location of the lowest point of the occultation
def Location(et, sv, when):
    Coords = np.ones(3)
    [tgopos, _] = spice.spkpos(sv.front, et-when, sv.fframe, 'NONE', sv.target)
    [mexpos, _] = spice.spkpos(sv.front, et-when, sv.fframe, 'NONE', sv.target)
    [sc2scvector,_] = spice.spkpos(sv.target, et-when, sv.fframe, 'NONE', sv.obs)
    displacement = np.linalg.norm(sc2scvector)
    sc2scunitvector = np.true_divide(sc2scvector, displacement)
    marsrad = spice.bodvrd(sv.front, 'RADII', 3)  # Extract the triaxial dimensions of Mars
    # For the ray that connects MEX and TGO, find the point on this ray that is closest to the Martian surface
    [nearestpoint,alt] = spice.npedln(marsrad[1][0], marsrad[1][1], marsrad[1][2], tgopos, sc2scunitvector)
    [radius, lon, lat] = spice.reclat(nearestpoint)
    lon = lon * (-180 / math.pi) # Rad -> Deg , frame inversion required (hence the negative 180)
    lat = lat * (-180 / math.pi)
    #Coords[0] = lon
    #Coords[1] = lat
    return lon,lat, displacement, nearestpoint, alt



# Normally SZA value is an unsigned integer, indicated the angle of whare the sun is. However this doesnt give
# any indication of what the local time would be (an angle could correspond to both AM and PM). Therefore the SZA is
# calculated twice to see if the angle is getting bigger or smaller at a later time. If smaller, then it is the AM (as
# the tangent point is moving towards the noon) Is idea of AM or PM is important for martian occultation as the
# ionosphere can vary wildly from sunrise-> sunset.
def SolarZenithAngles(et,nearestpoint, sv):

    subsolarpoint,_,_ = spice.subslr('INTERCEPT/ELLIPSOID',sv.front, et, sv.fframe, 'NONE', sv.target) # Where is the sun?
    sza = spice.vsep(subsolarpoint,nearestpoint) # angle between sun and sc
    latersubsolarpoint, _, _ = spice.subslr('INTERCEPT/ELLIPSOID', sv.front, et +30, sv.fframe, 'NONE', sv.target)# where is the sun later?
    latersza = spice.vsep(latersubsolarpoint, nearestpoint)

    if sza < latersza: #if the sun is moving away from the tangent point (PM)
        sza = sza * (180/math.pi)
        sza = sza*(-1) # negative SZA mean evening

    sza = sza * (180 / math.pi)# positive SZA means mornings

    return sza


def grazingangle(et,sv):

    profile,_, ht =  occgradient(sv.front, et, sv.fframe, sv.obs, sv.target) # find the shape of the occultation profile
    x,y,z = profile[0], profile[1], profile[2]
    lowestpoint = np.asarray([x[0],y[0],z[0]])
    highestpoint = np.asarray([x[ht-1],y[ht-1],z[ht-1]])
    maxsc2scvector = (highestpoint - lowestpoint) # find the vector from start to end of profile
    norm = np.linalg.norm(maxsc2scvector)
    endunitvector = maxsc2scvector/norm
    norm = np.linalg.norm(lowestpoint)
    startunitvector = lowestpoint/norm #provide unit vector for the lowestpoint, prepresenting the direction from Mars' center to the begining point
    angle = spice.vsep(startunitvector,endunitvector) #producing the anlge between Mars' center ->start & start-> end
    angle = angle * (180/math.pi)

    return angle


# Show the location of the occultation location on a map of mars
def charter(lon, lat,beg,stop,file_location):
    path_to_pic = file_location + '/images/2k_mars.jpg'
    raw_image = Image.open(path_to_pic)
    img = np.asarray(raw_image)# convert to array
    globe = ccrs.Globe(semimajor_axis=285000., semiminor_axis=229000., ellipse=None)
    crs = ccrs.PlateCarree(globe=globe)
    extent = (-895353.906273091, 895353.906273091, 447676.9531365455, -447676.9531365455) #adjustments of image for Robinson projection
    projection = ccrs.Robinson()
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=projection)
    ax.imshow(img, transform=crs, extent=extent)
    title = 'Occultation locations between '+beg[5:]+' and '+ stop[5:]
    plt.title(title)
    ax.plot(lon, lat, 'o',c='#bef9b9' , transform=ccrs.PlateCarree())
    plt.show()



# Produce a array containing all the coord (lon/lat), altitudes and SZAs for every meter along the Spacecraft-Spacecraft Vector
def producegeometrylamda(et,sv, when):
    [TGO, _] = spice.spkpos(sv.front, et-when, sv.fframe, 'NONE', sv.target)
    [MEX, _] = spice.spkpos(sv.front, et-when, sv.fframe, 'NONE', sv.obs)

    dist = math.floor(spice.vdist(TGO,MEX))
  
    # NEED TO PRODUCE VECTOR OF SZA AND HEIGHTS THAT ARE 'DIST' LONG [13242.9 m] *comp expensive 
    # start by making DIST length vector, 3 height, for every meter from mex to tgo
    # MAYBE FIND THE UNIT VECTOR AND ADD ONE IN ITS DIRECTION!!
    angleseparation = (spice.vsep(MEX, TGO))* (180/math.pi) # angle taken a mars center
    initialangle  = (spice.vsep(-MEX, (TGO-MEX))) * (180/math.pi)# angle taken at mars-MEX-tgo, that points to tgo. needed for the bending functions original starting angle  
    #script needs to work via periods of ray and not meters. [totalperiods is the main iterable, not meters]
    vacuumwavelength  = constants.c / 437.1e6
    scale = 1# scale =10, means we are itertating per 100 wavelenghts instead of 1000 (default 1000 because SPICE works in km)
    wavelengthsinameter = 1/vacuumwavelength
    a = wavelengthsinameter * dist * scale
    total1000periods = math.floor(a) # ~that many thousands of periods
    remainingdistance =(vacuumwavelength/scale)* ( (wavelengthsinameter * dist* scale) - total1000periods) # quanitfy the remaineder, this distance can
    # added later, this remaining portion is extreamly high altitude (near Target) and has no refractive effects. therfor simply added (km)
    #total1000periods = total1000periods.astype(int) 


    sc2sc = TGO - MEX
    norm = np.linalg.norm(sc2sc)
    unitsc2sc = sc2sc/(norm * vacuumwavelength * scale) #this needs to shrink if the repeatable expands
    points = np.empty([3,total1000periods])
    sza = np.empty([1,total1000periods])

    marsrad = spice.bodvrd(sv.front, 'RADII', 3)
    flatteningcoefficient = ( marsrad[1][0] - marsrad[1][2] ) / marsrad[1][0]
    equatorialradii = marsrad[1][0]
    # find direction of sun, it will not change much during the occultation. so only calc it once
    [SUN, _] = spice.spkpos(sv.front, et, sv.fframe, 'NONE', 'SUN')
    for i in range(total1000periods):
        point = MEX + (i * unitsc2sc) #move along ray, 1000 wavelength distance at a time (685 m). but unitsc2sc is in km...
        sza[0,i] = spice.vsep(SUN,point)
        points[:,i] = spice.recgeo(point, equatorialradii,flatteningcoefficient)
        points[0,i] = (points[0,i] * (-180 / math.pi))
        points[1,i] = (points[1,i] * (-180 / math.pi))
       
        
        #print((i/math.floor(total1000periods))*100)

    ray = np.concatenate((points,sza), axis=0)
  


    return ray, dist, total1000periods, remainingdistance

def producegeometrymeter(et,sv, when):
    [TGO, _] = spice.spkpos(sv.front, et-when, sv.fframe, 'NONE', sv.target)
    [MEX, _] = spice.spkpos(sv.front, et-when, sv.fframe, 'NONE', sv.obs)

    dist = math.floor(spice.vdist(TGO,MEX))
    print(dist)
    # NEED TO PRODUCE VECTOR OF SZA AND HEIGHTS THAT ARE 'DIST' LONG [13242.9 m] *comp expensive 
    # start by making DIST length vector, 3 height, for every meter from mex to tgo
    # MAYBE FIND THE UNIT VECTOR AND ADD ONE IN ITS DIRECTION!!
    angleseparation = (spice.vsep(MEX, TGO)) # angle taken a mars center
    initialangle  = (spice.vsep(-MEX, (TGO-MEX))) * (180/math.pi)# angle taken at mars-MEX-tgo, that points to tgo. needed for the bending functions original starting angle  
    #script needs to work via periods of ray and not meters. [totalperiods is the main iterable, not meters]
    
    scale = 0.5 # scale =10, means we are itertating per 100 wavelenghts instead of 1000 (default 1000 because SPICE works in km)
    
    
    dist = math.floor(dist) # km

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
       
        
        print((i/math.floor(dist))*100)

    ray = np.concatenate((points,sza), axis=0)

    #plt.plot(angleprogression[0,:], ray[2,:])
    #plt.show()

    # ray is in lat/lon/alt + sza and xyzpoints is cartesian, both describe the same thing
    return ray, dist , unitsc2sc,angleseparation, initialangle, MEX,TGO,  xyzpoints, angleprogression



def flatbending(xyzpoints,initialangle, sv, MEX,TGO):
    #form a coordinate system where tgo is @ y=0 and x= (5000 +norm), Mar's Barrycenter being @ [5000,0]

    #initialise non-global variables
    miniray = np.zeros(10)
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
    for i in range(length):
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
    direction = unitmex.dot(rotationvector) * unit #make a 2d vector coming from MEX YOU DO NOT KNOW WHAT WAY THIS IS ROTATING


    iterationcount =0
    #while iterationcount<100:
    print( "Finding Bending Angle (", str(iterationcount) ,"% Complete)")

    #some function based on miss (miss begins with )

    stage =0
    t=0 # index counter for steps along the ray
    while stage < 2: #==0first unit, so move two units. ==1 propergate step by step. ==2 exit and analyse entire path 
        #take the alt at 10 possitions across this unit
        if stage==0:
            for k in range(10):
                point = mex + (k *(direction/10))
                #_,_,miniray[k] = spice.recgeo(point, equatorialradii,flatteningcoefficient)
                miniray[k] =  np.linalg.norm(point) - 3389 #average radii of mars
            N0 = findrefractivity(miniray,10)
            raystep[:,t] = point #save the last location
            t=t+1
            stage =stage+1

        if stage==1:
            for k in range(10):
                point = raystep[:,t-1] + (k *(direction/10))
                #_,_,miniray[k] = spice.recgeo(point, equatorialradii,flatteningcoefficient)  #THIS ONLY WORKS IN 3D
                # IMPLEMENTING MARS AS A SIMPLE CIRCLE OF AVERAGE 3389 KM RADIUS, !THIS WILL BE UPDATED TO ELLIPSE!
                miniray[k] =  np.linalg.norm(point) - 3389
            raystep[:,t] = point#9 is the end of the unit, and refraction always happens relative to the center of refractivity, so rotate off this vector
            N1 = findrefractivity(miniray,10)

            if point[1] < 0: #if the position drops below the x axis
                stage = stage+1 #increase the stage value so the while loop is exited

            # if N0 !=1:
            #     print('querycode')

            currenty = raystep[1,t]
            print('Current Y possition is', currenty)
            r= N0/N1 #NEED MORE PRECISION
            if r==1: #only bend when there is a refractive gradient between consecutive air volumes[NO CHANGE IN DIRECTION]
                t=t+1
                N0=N1
                continue

            #find the angle between the unit (air volume) boarder and the current direction
            unitrotationaxis = -(raystep[:,t]/np.linalg.norm(raystep[:,t]))
            #unitdirection = direction #MAYBE ALTERING UNITS WILL EFFECT THIS * OR / BY UNIT, CANT FIGURE OUT NOW, OK WHEN UNIT =1
            DotProduct=  np.dot(unitrotationaxis,direction)
            AngleofIncidence = (math.pi/2)- np.arccos(DotProduct) #angle it enters the next air volume
            #simple snell law to find the bending angle (should be tiny angle)     
            AngleofRefraction = np.arcsin(r * np.sin(AngleofIncidence))
            # THIS IS NOT EXACTLY WHAT THE TURN IN DIRECTION IS, NEED TO THINK ABOUT
            rotateby = 300*((AngleofIncidence - AngleofRefraction))#+ve =clockwise, -ve=anticlockwise

            INCIDENCEDEGREES = np.degrees(AngleofIncidence)
            REFRACTIONDEGREES = np.degrees(AngleofRefraction)
            ROTATIONDEGREES = np.degrees(rotateby)

            #an if statement is required, if this 
            if ROTATIONDEGREES>1 or ROTATIONDEGREES<-1:
                print('stophere, u r bending to much')
            rotationvector = np.array(( (np.cos(rotateby), -np.sin(rotateby)),
            (np.sin(rotateby),  np.cos(rotateby)) ))
            direction = direction.dot(rotationvector)

            #if np.linalg.norm(direction)< 0.08:# assuming unit is set to 0.1
                # print(' your movement direction is shrinking')

            t=t+1
            N0=N1

    error = np.zeros(t)   
    g=0

    #error = finderror(raystep, UnrefractedRay, initialtheta)

    miss = error[-1] # +ve miss needs a clockwise rotation, 
    iterationcount=iterationcount+1



    fig, ax = plt.subplots()
    plt.plot(barry[0],barry[1],'x' )
    plt.annotate('$\u2642$', (barry[0],barry[1]), fontsize = 20)
    marsradii = plt.Circle((0,0),3389,color = 'red', fill=False)
    ax.add_artist(marsradii)
    plt.plot(mex[0],mex[1], '.')
    plt.annotate('MEX',(mex[0],mex[1]))
    plt.plot(tgo[0], tgo[1], 'o')
    plt.annotate('TGO',(tgo[0], tgo[1]) )
    plt.plot(UnrefractedRay[0],UnrefractedRay[1], ':')
    plt.annotate('Distance = %ikm'%UnrefractedDistance, (UnrefractedRay[0,length//2],UnrefractedRay[1,length//2]), fontsize = 8)
    plt.plot(raystep[0], raystep[1], ':')
    plt.gca().set_aspect('equal', adjustable='box')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    plt.show()

    print('stophere')


def finderror(refract, unrefract, angle):
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

    shortenrefract = refract[0,:]
    shortenrefract = shortenrefract[shortenrefract !=0]
    error = np.zeros(len(shortenrefract)-1) #initailise a 2d vector

    for i in range(len(shortenrefract)-1):
        ox, oy = unrefract[:,0] #extract MEX as origin
        px, py = refract[:,i]
        rotx = ox + math.cos(-newangle) * (px - ox) - math.sin(-newangle) * (py - oy)#rotate each point by initial starting angle
        error[i] = rotx

    error = error - unrefract[0,0] #subtract the origin x-value from every value in error

    plt.plot(error)
    plt.ylabel('Error from Unrefracted Path (km)')
    plt.xlabel('Propergation Along Ray Path MEX->TGO (km)')
    plt.show()
    return error

# Find the average refractivty of volume descriped by ray
def findrefractivity(alt, step): 
    #must convert these to 32decimal
    ionoresidual = atmosphere.iono(alt,step)
    neutralresidual = atmosphere.neutral(alt,step)
    n = np.sum(1 + (ionoresidual + neutralresidual))/step #produce the integral refractive index and divide by lenght to find average
    return n



#the issue must be here then, the function above makes sense.
def doppler(residual, totalperiods, dist, remaining): # the results of this is changing for the the number of iterations
    transmitfrequency = 437.1e6
    #dist needs to be altered, how many wavelenghts fit into the whole dist, that needs to be the iterable. wavelength is 61 cm. this is too small
    scale = 1 # scale =10, means we are itertating per 100 wavelenghts instead of 1000
    vacuumwavelength  = (constants.c / transmitfrequency) /scale
    total1000periods = ((1/vacuumwavelength) * dist * scale) # to test if this is the rounding error in the geometricdistasnce 
    wavelength = np.ones([1,totalperiods])

    for i in range(totalperiods):
        wavelength[0,i] = (vacuumwavelength) * residual[0,i] #resisdual is only dist long

    electricdistance = np.sum(wavelength)
    geometricdistance = (vacuumwavelength*total1000periods)/scale # now this has not got a rounding error
    delta = (electricdistance - geometricdistance + remaining) 
    
    dopplershift = ((delta/geometricdistance) * transmitfrequency) 
    return electricdistance , geometricdistance, dopplershift