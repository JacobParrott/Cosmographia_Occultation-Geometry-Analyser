import numpy as np
import spiceypy as spice
from scipy import constants
from scipy.interpolate import InterpolatedUnivariateSpline
import math
import array as arr
import shapely
#from decimal import *
from tqdm import tqdm
import time as timer
from mpmath import *
mp.dps = 60 
#getcontext().prec = 35 #set the precision for the bending in function "flatbending" 



import JSONtemplate as JS
import pandas as pd
import json
from PIL import Image
import matplotlib.pyplot as plt
import cartopy.crs as ccrs #import the coordinate refernece system
from scipy.spatial.transform import Rotation as R

import atmosphere_highprecision


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
    for i in tqdm(range(dist)):
        xyzpoint = MEX + (i * unitsc2sc) #move along ray, 1000 wavelength distance at a time (685 m). but unitsc2sc is in km...
        xyzpoints[:,i] = xyzpoint
        sza[0,i] = spice.vsep(SUN,xyzpoint)
        angleprogression[0,i] = (spice.vsep( xyzpoint, MEX)) * (180 / math.pi)
        points[:,i] = spice.recgeo(xyzpoint, equatorialradii,flatteningcoefficient)
        points[0,i] = (points[0,i] * (-180 / math.pi))
        points[1,i] = (points[1,i] * (-180 / math.pi))
       
        
        

    ray = np.concatenate((points,sza), axis=0)

    #plt.plot(angleprogression[0,:], ray[2,:])
    #plt.show()

    # ray is in lat/lon/alt + sza and xyzpoints is cartesian, both describe the same thing
    return initialangle, MEX,TGO,  xyzpoints



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
    UnrefractedRay = np.zeros((2,length))
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

    unit = 0.5 # in km
    
    
    rotationvector = np.array(( (np.cos(initialtheta), -np.sin(initialtheta)),
                                (np.sin(initialtheta),  np.cos(initialtheta)) ))

    #get unit vecotr of -MEX (then add this vecotr to MEX for each alt calcultation)
    unitmex = -mex/barry2mex #unit direction (2d)
    initialdirection = unitmex.dot(rotationvector) * unit #make a 2d vector coming from MEX YOU DO NOT KNOW WHAT WAY THIS IS ROTATING


    iterationcount =0
    #while iterationcount<100:
    #print( "Finding Bending Angle (", str(iterationcount) ,"% Complete)")
    errorstore = np.zeros((11,5000000))
    miss= inf
    missangle= inf 
    
    while iterationcount <100:
        if iterationcount ==0:
            direction = initialdirection
        else:
            # the initail direction must be rotated by a 10th of the miss at the end
            missangle = missangle/1
            missrotationvector = np.array(( (mp.cos(missangle),- mp.sin(missangle)),
                                            (mp.sin(missangle), mp.cos(missangle)) ))
            direction = initialdirection.dot(missrotationvector)
           

            #to make the produced direction reach mpmath dp level
            # a= mp.cos(missangle) ; b = -mp.sin(missangle) #ORIGINAL DID NOT GAVE THIS NEG, SO IT HAS BEEN MADE POSITIVE
            # c = mp.sin(missangle) ; d = mp.cos(missangle)

            # rotationvector = np.array([[a,b],[c,d]])
            # direction = direction.dot(rotationvector)

            initialdirection = direction 
        Nstore = np.zeros(2000000)
        turningcounter=0
        progress= [0,0]
        stage =0
        t=0
        tic  = timer.perf_counter()
        
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
                    raystep[0,t] = point[0] ; raystep[1,t] = point[1] #save the last location
                    t=t+1
                    stage =stage+1

                if stage==1:
                    for k in range(subgroupsize):
                        subvalue = ((k+1) *(direction/subgroupsize))
                        point = raystep[:,t-1] + subvalue #am i double counting the end of the last and the start of the next?
                        #_,_,miniray[k] = spice.recgeo(point, equatorialradii,flatteningcoefficient)  #THIS ONLY WORKS IN 3D
                        # IMPLEMENTING MARS AS A SIMPLE CIRCLE OF AVERAGE 3389 KM RADIUS, !THIS WILL BE UPDATED TO ELLIPSE!
                        miniray[k] =  np.linalg.norm(point) - 3389 #provides the alt, this can have less precision
                    raystep[0,t] = point[0] ; raystep[1,t] = point[1] #save the last location
                    N1 = findrefractivity(miniray,subgroupsize)

                    if point[1] < 0: #if the position drops below the x axis
                        stage = stage+1 #increase the stage value so the while loop is exited


                    currenty = raystep[1,t]
                    progress[1] = mex[1] - currenty #this value will be increasing from 0 -> Mex height
                    increment = np.floor(progress[1])- np.floor(progress[0]) #only when 
                    if increment ==1 :
                        pbar.update(1)
                    progress[0] = progress[1]

                    if abs(N0) < 1e-50:#catch for precision errors e-50 is basically the entire propagation
                        Nstore[t] = N1
                        t=t+1
                        N0=N1
                        continue
                    
                    #ensure all the turns are being implemented

                    #print('Current Y possition is', currenty) #this is alt, so is ~3389 km smaller than the vector
                
                    numorator = mpf(N0)+mpf(1)
                    denominator = mpf(N1)+mpf(1)
                    rbending = mp.fdiv(numorator,denominator) #average would just add 1 to total N
                    if t==5000: #only bend when there is a refractive gradient between consecutive air volumes[NO CHANGE IN DIRECTION]
                        t=t+1
                        N0=N1
                        continue

                    #this section is only reached if a turning is going to happen
                    # ,testing to see if there is 10 X less turing if the units are 10X smaller -TRUE
                    # TEST = float(rbending)
                    # if TEST != 1:
                    #     turningcounter = turningcounter+1


                    # !! NOW WITH PRECISION !!

                    #find the angle between the unit (air volume) boarder and the current direction
                    unitrotationaxis = raystep[:,t]/np.linalg.norm(raystep[:,t])
                    #unitdirection = direction #MAYBE ALTERING UNITS WILL EFFECT THIS * OR / BY UNIT, CANT FIGURE OUT NOW, OK WHEN UNIT =1
                    DotProduct=  fdot(unitrotationaxis,direction)
                    AngleofIncidence = (math.pi/2)- mp.acos(DotProduct) #angle it enters the next air volume
                    #simple snell law to find the bending angle (should be tiny angle)     
                    AngleofRefraction = mp.asin(rbending * mp.sin(AngleofIncidence))
                    # THIS IS NOT EXACTLY WHAT THE TURN IN DIRECTION IS, NEED TO THINK ABOUT
                    rotateby = ((AngleofIncidence - AngleofRefraction))#+ve =clockwise, -ve=anticlockwise

                    INCIDENCEDEGREES = mp.degrees(AngleofIncidence)
                    REFRACTIONDEGREES = mp.degrees(AngleofRefraction)
                    ROTATIONDEGREES = mp.degrees(rotateby)

                    predirection = direction

                    #an if statement is required, if this 
                    if ROTATIONDEGREES>1 or ROTATIONDEGREES<-1:
                        print('stophere, u r bending to much')
                    a= mp.cos(rotateby) ; b = -mp.sin(rotateby)
                    c = mp.sin(rotateby) ; d = mp.cos(rotateby)

                    rotationvector = np.array([[a,b],[c,d]])
                    direction = direction.dot(rotationvector)



                    # rotationvector = matrix([[ mp.cos(rotateby), -mp.sin(rotateby)],
                    #                           [mp.sin(rotateby),  mp.cos(rotateby)]])
                    # direction = fdot(direction,rotationvector)

                    if direction[1] != predirection[1] :
                        #no bending has occured due to precision error
                        #if the last direction is different to the current one then a turning has occured.
                        #do we move forward with the same precision
                        turningcounter = turningcounter+1

                    N0=N1

                    #store N1 to calc electric distance 
                    Nstore[t] = N1
                    t=t+1
        
        unit_initial = initialdirection/ np.linalg.norm(initialdirection)
        dot_product = np.dot(unit_initial,unitmex)
        FinalBendingAngle = mp.acos(dot_product)

        #from the refractive profile from MEX ->TGO, calc the intergral of the change in wavelength to aquire doppler
        Nstore =  Nstore[Nstore !=0]
        Nstore =  1-Nstore #convert N deltas into values of N

       
     
       

        error = np.zeros(t)   
        #print("Number of turns:", turningcounter)
        error = finderror(raystep, UnrefractedRay, initialtheta, iterationcount)
        miss = error[-1] # +ve miss needs a clockwise rotation

        errorstore[iterationcount,:t-1] = error

        missangle = mp.asin(miss/ UnrefractedDistance)
        toc  = timer.perf_counter()
        passingtime = toc-tic
        print(' miss =', format(miss*1000, '.5f') ,'m || Angle =', nstr(mp.degrees(FinalBendingAngle + initialtheta),5) ,
                            '° || Speed =',passingtime,' Sec \n', sep = " ", end= " ", flush =True)
        
        #print('ANGLE miss =', mp.degrees(missangle) ,'Degrees')
        if abs(miss) <1e-5: #is the miss smaller than 10 cm?
            #ploterrortraces(errorstore,t)
            print("stophere")
            break
        iterationcount=iterationcount+1

    #Expensive plotting (7sec)
                # fig, ax = plt.subplots()
                # plt.plot(barry[0],barry[1],'x' )
                # plt.annotate('$\u2642$', (barry[0],barry[1]), fontsize = 20)
                # marsradii = plt.Circle((0,0),3389,color = 'red', fill=False)
                # ionoradii = plt.Circle((0,0),3389+130,color = 'blue', fill=False)
                # ax.add_artist(marsradii)
                # ax.add_artist(ionoradii)
                # plt.plot(mex[0],mex[1], '.')
                # plt.annotate('MEX',(mex[0],mex[1]))
                # plt.plot(tgo[0], tgo[1], 'o')
                # plt.annotate('TGO',(tgo[0], tgo[1]) )
                # plt.plot(UnrefractedRay[0],UnrefractedRay[1], ':')
                # plt.annotate('Distance = %ikm'%UnrefractedDistance, (UnrefractedRay[0,length//2],UnrefractedRay[1,length//2]), fontsize = 8)
                # plt.plot(raystep[0], raystep[1], ':')
                # plt.gca().set_aspect('equal', adjustable='box')
                # ax.set_xlabel('X')
                # ax.set_ylabel('Y')
                # plt.plot(raystep[0,0::500],raystep[1,0::500], 'x', markersize = 12)
                # plt.show()


    ElectricDistance = np.sum(Nstore) #N * wavelengths in a km (UNITS MUST BE KEPT THE SAME AS STEP-SIZE)

    return FinalBendingAngle, ElectricDistance
    print('stophere')


def finderror(refract, unrefract, angle,ittcount):
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

    for i in range(len(shortenrefract)-1): #Rotate the propagation so it heads vertically downwards
        ox, oy = unrefract[:,0] #extract MEX as origin
        px, py = refract[:,i]
        rotx = ox + math.cos(-newangle) * (px - ox) - math.sin(-newangle) * (py - oy)#rotate each point by initial starting angle, so it goes downwards
        error[i] = rotx

    error = error - unrefract[0,0] #subtract the origin x-value from every value in error
    plottingerror = error * -1000 #convert to m

    ########
    #interpolate the error so the first and second derivatives can be found
    # x = np.array(range(len(shortenrefract)-1))
    # smoothresults = InterpolatedUnivariateSpline(x, plottingerror)
    # prime = smoothresults.derivative(n=1)
    # doubleprime = smoothresults.derivative(n=2)
    # finerx = np.linspace(0,len(x),len(x))
    # firstderivative= prime(finerx)
    # secondderivative= abs(doubleprime(finerx))

    # fig, ax = plt.subplots(3)
    
    # ax[0].plot(plottingerror, 'r')
    # ax[1].plot(firstderivative)
    # ax[2].plot(secondderivative,'g')
    # ax[0].set_title('Error from Unrefracted Path', loc = 'left')
    # ax[0].set_ylabel('Posistion from Unrefracted Path (m)')
    # ax[1].set_title('Error`', loc = 'left')
    # ax[1].set_ylabel('Gradient from Unrefracted Path')
    # ax[2].set_title('Error``', loc = 'left')
    # ax[2].set_ylabel('Cange Gradient from Unrefracted Path')
    # ax[0].set_xlabel('Propergation Along Ray Path MEX->TGO (km)')
    # plt.show()
    ########

    # just want to see the error plot (not differentiate a huge dataset)
   
    plt.plot(plottingerror, 'r')
    plt.title('Error from Unrefracted Path', loc = 'left')
    plt.ylabel('Posistion from Unrefracted Path (m)')
    plt.xlabel('Propergation Along Ray Path MEX->TGO (km)')
    plt.show()

    return error

# Find the average refractivty of volume descriped by ray
def findrefractivity(alt, step): 
    #must convert these to 32decimal
    ionoresidual = atmosphere_highprecision.iono(alt,step)
    neutralresidual = atmosphere_highprecision.neutral(alt,step)
    sumiono = fsum(ionoresidual) ; sumneut = fsum(neutralresidual)
    sumtotal = fadd(sumiono,sumneut)
    n = fdiv(sumtotal,step) #produce the integral refractive index and divide by lenght to find average
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

def ploterrortraces(errorstore,t):

    fig, ((ax1,ax2),
        (ax3,ax4),
        (ax5,ax6),
        (ax7,ax8),
        (ax9,ax10)) = plt.subplots(5,2)

    fig.suptitle("Error from Unrefracted Path by Iteration (m)")
    ax1.plot(-1000*errorstore[0,:t-1]) ; ax1.axhline(0, c='r', ls='--', lw=.5); ax1.annotate("1)",(0,0))
    ax2.plot(-1000*errorstore[1,:t-1]) ; ax2.axhline(0, c='r', ls='--', lw=.5); ax2.annotate("2)",(0,0))
    ax3.plot(-1000*errorstore[2,:t-1]) ; ax3.axhline(0, c='r', ls='--', lw=.5); ax3.annotate("3)",(0,0))
    ax4.plot(-1000*errorstore[3,:t-1]) ; ax4.axhline(0, c='r', ls='--', lw=.5); ax4.annotate("4)",(0,0))
    ax5.plot(-1000*errorstore[4,:t-1]) ; ax5.axhline(0, c='r', ls='--', lw=.5); ax5.annotate("5)",(0,0))
    ax6.plot(-1000*errorstore[5,:t-1]) ; ax6.axhline(0, c='r', ls='--', lw=.5); ax6.annotate("6)",(0,0))
    ax7.plot(-1000*errorstore[6,:t-1]) ; ax7.axhline(0, c='r', ls='--', lw=.5); ax7.annotate("7)",(0,0))
    ax8.plot(-1000*errorstore[7,:t-1]) ; ax8.axhline(0, c='r', ls='--', lw=.5); ax8.annotate("8)",(0,0))
    ax9.plot(-1000*errorstore[8,:t-1]) ; ax9.axhline(0, c='r', ls='--', lw=.5); ax9.annotate("9)",(0,0))
    ax10.plot(-1000*errorstore[9,:t-1]) ; ax10.axhline(0, c='r', ls='--', lw=.5); ax10.annotate("10)",(0,0))
    
    plt.show()