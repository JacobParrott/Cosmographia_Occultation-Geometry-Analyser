import numpy as np
import spiceypy as spice
from scipy import constants
from scipy.interpolate import interpolate
import math
import array as arr
import shapely

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
    print(dist)
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
       
        
        print((i/math.floor(total1000periods))*100)

    ray = np.concatenate((points,sza), axis=0)
    print('stop here')


    return ray, dist ,angleseparation, initialangle, total1000periods, vacuumwavelength, remainingdistance

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
    
    scale = 0.1 # scale =10, means we are itertating per 100 wavelenghts instead of 1000 (default 1000 because SPICE works in km)
    
    
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



def cartesianbending(xyzpoints,initialangle, sv, MEX,TGO):
    #form a coordinate system where tgo is @ y=0 and x= (5000 +norm), Mar's Barrycenter being @ [5000,0]

    #initialise non-global variables
    miniray = np.zeros(10)
    raystep = np.zeros((2,10000))# create a large array to populate and then shrink later

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


    #psuedo
    #.decide what the unit is (probably iterable),start propergation with the initial angle in one unit direction (doesnt have to be 1 km),
    #. feed each unit/10 into an alt finder and calc the avg N for this unit
    #. simple snell to find new angle with n0/n1 (does it need entry angle?)
    #.find barry2ray possition, and the bend will be from the normal of that (simple)
    #.rotate current vector up(iono) or down(neutral) minutely
    #.save coords at each unit step
    # . propergate
    #.exit when y = ~0 (bare in mind the nuetral and really fuck up at low alts), this should show in a plot anyways
    #. check variance with an altered ploterrorvariation() function
    #.calc miss, x value will usually be greater than tgo(x) because it will mostly bend upwards

    initialtheta = -(spice.vsep(MEX-TGO, MEX)) #this will produce and angle that is likly not going to be exactly on 
    #the original propagation path, you compare to this if there is a drifting error, as both this and the resultant refracted ray 
    # have the same bias error. THIS ANGLE IS BENDING ANTICLOCKWISE IN THIS FRAME (BENDING UPWARDS)
    nicetohave = np.degrees(initialtheta)
    unit = 1 # in km
    rotationvector = np.array(( (np.cos(initialtheta), -np.sin(initialtheta)),
               (np.sin(initialtheta),  np.cos(initialtheta)) ))

    #get unit vecotr of -MEX (then add this vecotr to MEX for each alt calcultation)
    unitmex = -mex/barry2mex #unit direction (2d)
    direction = unitmex.dot(rotationvector) * unit #make a 2d vector coming from MEX YOU DO NOT KNOW WHAT WAY THIS IS ROTATING

    stage =0
    t=0 # index counter for steps along the ray
    while stage < 2: #==0first unit, so move two units. ==1 propergate step by step. ==2 exit and analyse entire path 
        #take the alt at 10 possitions across this unit
        if stage==0:
            for k in range(10):
                point = mex + (k *(direction/10))
                #_,_,miniray[k] = spice.recgeo(point, equatorialradii,flatteningcoefficient)
                miniray[k] =  np.linalg.norm(point) - 3389
            N0 = findrefractivity(miniray,10)
            raystep[:,t] = point #save the last location
            t=t+1
            stage =stage+1

        if t==5440:
            print('stophere')
        
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

            #DEBUGGING CATCHES: 
            r= N0/N1
            if r==1: #only bend when there is a refractive gradient between consecutive air volumes
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

            t=t+1
            N0=N1

    error = np.zeros(t)   
    g=0

    #INSERT THE SOLUTION HERE:
    #method4
    error = finderror(raystep, UnrefractedRay, initialtheta)

    #smoothraystep  = interpolate.interp1d(raystep[0],raystep[1], kind = 'linear')
    #smoothUnrefractedRay  = interpolate.interp1d(UnrefractedRay[0],UnrefractedRay[1], kind = 'linear') 
    # i=0
    # raysteparea = np.trapz(raystep[0:6482], axis =0)
    # unrefractedarea = np.trapz(UnrefractedRay[0:6482], axis =0)
    # #for i in range(6482):

    # error[i] = raysteppoint[0]-UnrefractedRaypoint[0] + 0.25
    #x_progress  = np.floor(range(int(np.floor(max(raystep[0,:])) - np.floor(min(raystep[0,:])))) + raystep[0,0])
    # raystep_xy = np.zeros((2,np.size(x_progress)))
    #UnrefractedRay_xy = np.zeros((2,np.size(x_progress)))
    error = np.zeros(t) 
    # smallerraystep_x = raystep[0,:]
    # smallerraystep_x = smallerraystep_x[smallerraystep_x != 0]
    # roundedraysteps_x = np.floor(smallerraystep_x)
    # roundedunrefracted_x = np.floor(UnrefractedRay[0,:])

    # #rounded_x = np.floor(UnrefractedRay[0,:])
    # #analyse at each raystep to see the error, you can do an x value match and then correct for the delay
    #for i in range(t):

        # raystep_x = raystep[0,:]
        # UnrefractedRay_x = UnrefractedRay[0,:]
        # idxray = np.searchsorted(UnrefractedRay_x , raystep_x[i], side ='left')
        # if idxray > UnrefractedRay_x[-1]:
        #     break
        # UnrefractedRay_xy = UnrefractedRay[:,idxray]
        # #idxUnrefracted = np.searchsorted(roundedunrefracted_x,  x_progress[i], side ='left') # and raystep_x[0,:] > x_progress[i] ]# to deal with the high precision numbers
        # if idxray ==[]:
        #     break
        
        #ray_y = raystep[1,idxray]
        #unrefracted_y = UnrefractedRay[1,idxUnrefracted]
        #error[i] = ray_y-unrefracted_y

        # index_xUnrefracted= np.where(UnrefractedRay[0,:] < (x_progress[i] + 1) and UnrefractedRay[0,:] > x_progress[i] ) # it is not found because Ray has many decimal points
        # index_xRefracted = np.where(raystep[0,:] < (x_progress[i] +1) and raystep[0,:] > x_progress[i]) #maybe add a -+ limiter aound it?
        # indexU = index_xUnrefracted[1][0]
        # UnrefractedRaypoint = UnrefractedRay[:,indexU]
        # indexR = index_xRefracted[1][0]
        # raysteppoint = raystep[:,indexR]
        


    #     error[i] = raysteppoint[0]-UnrefractedRaypoint[0] + 0.25
    # #     #PossiblePoints = PossiblePoints[PossiblePoints > x_progress[i] ]# splitting the conditional across two lines due to a ambiguity error
    # #     #PossiblePoints1 = PossiblePoints[PossiblePoints !=0]

    # error = error[error !=0]

    # #LETS SMOOTH THE ERRORS TO FIX YOUR ROUNDING ERROR
    # window_size = 10
    # i = 0
    # moving_averages = []
    # while i < len(error) - window_size + 1:
    #     this_window = error[i : i + window_size]
    #     window_average = sum(this_window) / window_size
    #     moving_averages.append(window_average)
    #     i += 1

    # plt.plot( moving_averages , '.')
    # plt.ylabel('Error')
    # plt.xlabel('X progression')
    # plt.show()



    #     # #FUCK IT EVEN OLDER VERSION
    #     # #FIND VALUE IN RAYSTEP THAT SUITS THE XPROGRESS CONDITIONS
    #     # raystep_x = raystep[0,:]
    #     # PossiblePoints = raystep_x[raystep_x < (x_progress[i]+1)]# and raystep_x[0,:] > x_progress[i] ]# to deal with the high precision numbers
    #     # PossiblePoints = PossiblePoints[PossiblePoints > x_progress[i] ]# splitting the conditional across two lines due to a ambiguity error
    #     # PossiblePoints1 = PossiblePoints[PossiblePoints !=0]
    #     # if PossiblePoints.size ==0:#very rarely there maybe be no value for this value of x
    #     #     continue #we might have got to the end 
    #     # SearchFor = PossiblePoints1[0]# only take the first value
    #     # index = np.where(raystep_x == SearchFor)
    #     # Index1 = index[0][0]
    #     # raystep_xy[:,i] = raystep[:,Index1]

    #     # #FIND VALUE IN UNREFRACTED THAT SUITS THE XPROGRESS CONDITIONS
    #     # UnrefractedRay_x = UnrefractedRay[0,:]
    #     # PossiblePoints = UnrefractedRay_x[UnrefractedRay_x< (x_progress[i]+1)]# to deal with the high precision numbers
    #     # PossiblePoints = PossiblePoints[PossiblePoints > x_progress[i] ]# splitting the conditional across two lines due to a ambiguity error
    #     # PossiblePoints2 = PossiblePoints[PossiblePoints !=0]
    #     # if PossiblePoints.size ==0:
    #     #     continue
    #     # SearchFor = PossiblePoints2[0]# only take the first value
    #     # index = np.where(UnrefractedRay_x == SearchFor)
    #     # Index2 = index[0][0]
    #     # UnrefractedRay_xy[:,i] = UnrefractedRay[:,Index2]
        
    #     # # OLD METHOD

    #     #error[i] = np.linalg.norm(raystep_xy[:,i]-UnrefractedRay_xy[:,i]) #probs alter the format of 'index'
    #     g=g+1 #only do this if it set!!

    # #error = error[error < 0.1 ]#take every 10th value]
    # #x_progress = x_progress[0::10]
    # plt.plot( error , '.')
    # plt.ylabel('Error')
    # plt.xlabel('X progression')
    # plt.show()
    
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
    plt.plot(raystep[0], raystep[1], '.')
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
        #roty = oy + math.sin(angle) * (px - ox) + math.cos(angle) * (py - oy)

        error[i] = rotx

    error = error - unrefract[0,0] #subtract the origin x-value from every value in error

    plt.plot(error)
    plt.show()
    return error

# Find the average refractivty of volume descriped by ray
def findrefractivity(alt, step): 
    ionoresidual = atmosphere.iono(alt,step)
    neutralresidual = atmosphere.neutral(alt,step)
    n = np.sum(1 + (ionoresidual + neutralresidual))/step #produce the integral refractive index and divide by lenght to find average
    return n



#iterate coarsly through the angle of separation until TGO is met, both dl and angle change. 
#spacial varies with change in theta (2pitheta)
def bending( ray, initialdirection, MEX,  angleseparation, angleprogression, initialangle, xyzpoints,sv):
    alt = ray[2,:]
   
    #THIS IS AN ATMOSPHERE TEST
    # ionoresidual = atmosphere.iono(ray,np.size(ray,1))
    # neutralresidual = atmosphere.neutral(ray,np.size(ray,1))
    # residual = 1 + (ionoresidual + neutralresidual)

    # plt.plot(angleprogression[0,:],residual[0,:])
    # plt.title("Refractive Index through Propergation of Ray")
    # plt.xlabel("MEX -> Possition angle (degrees)")
    # plt.ylabel("Refractive Index")
    # plt.show()

    rawvector = xyzpoints[:,np.size(xyzpoints,1)-1]-xyzpoints[:,0]
    norm = np.linalg.norm(rawvector)
    rayunit = rawvector/norm # THIS IS SLIGHTLY DIFFERENT TO THE RAY, THERE ARE PYTHON ROUNDING ERRORS



    marsrad = spice.bodvrd(sv.front, 'RADII', 3)
    flatteningcoefficient = ( marsrad[1][0] - marsrad[1][2] ) / marsrad[1][0]
    equatorialradii = marsrad[1][0]
    direction = rayunit
    # everything must have a alt/ equiv dtheta
    poss =np.zeros([3,100000])
    straightray = np.zeros([3,100000])
    poss[:,0] = MEX
    straightray[:,0] = MEX
    iterationmax = 0
    miss=0
    glitchcount = 0
    totalarcdistance = math.floor(alt[1] * angleseparation) # in km
    angleseparation = angleseparation * (180 / math.pi) #we want in degrees now
    #NOTE TO SELF, AT FIRST THE DIFFERNCE IN N MIGHT BE TO SMALL IN THE HIGH ALT. MAYBE JUMP TO WHERE THE IONO BEGINS THEN START THIS LOOP

    while iterationmax <=100:
        dtheta = np.zeros(100000)# must be re-initialised as the vectors are filtered for 0 towards the end of this loop
        possalt =np.zeros(100000)
        straightpossalt =np.zeros(100000)
        #alter initial angle slightly according to miss
        #if miss ~0 dont alter angle, if miss +ve then bend down, (you will likely be bending up initially due to tropo)

        if miss ==0:
            print('.')# continue will make a 10-line loop with while
        else:
            angle = np.radians(miss * 0.03) # a miss of one km (0.3km = 2 degrees )
        #alter the initail starting angle
            unitposs =    MEX / np.linalg.norm(MEX) #vector going upwards radially if n0>n1
            unitdirection = initialdirection/ np.linalg.norm(initialdirection)
            rotationaxis = np.cross(unitposs, unitdirection) #if direction is going rightwards, this rotation vector is going away from the viewer
            rotationvector = rotationaxis * angle  #must be in degrees
            rotation = R.from_rotvec(rotationvector)
            direction = rotation.apply(initialdirection) # produce the slightly re-angled new direction (tiny difference)



        step = 10  # this will scale down to 0.001 (1m)
        interval = 1
        miniray = np.zeros([3,10])
        constantdirection = direction
        for p in range(100000//step): #10,000 is the max distance, this will never be reached
            if p == 0: #You need a before n and after n so for the begining go forwards two to get a n0 and n1
                point = poss[:,0] +0 # regeo needs non-strided arrays, this achieved by adding 0 (cheat)
                _,_, possalt[p] = spice.recgeo(point, equatorialradii,flatteningcoefficient)
                dtheta[p] = (spice.vsep( point, MEX))* (180 / math.pi) #should give 0 degrees
            
                for i in range(step):
                    point = poss[:,0] + (interval * i * direction) #make a small 'step' long 3d vector
                    miniray[:,i] = spice.recgeo(point, equatorialradii,flatteningcoefficient) #(lon; lat; alt)
                n0 = findrefractivity(miniray,step)
                
            elif p==477:
                print('inspect from here')    
                
            else:
                n0=n1

            for i in range(step):
                point = poss[:,p] + (interval * i * direction) #make a small 'step' long 3d vector
                miniray[:,i] = spice.recgeo(point, equatorialradii,flatteningcoefficient)
            n1 = findrefractivity(miniray,step)
            

            #VECTOR CLACS [dont understand inuitivly, boiler plate] maybe this breaks if no r
            #SHOULD HIT THE IONO AT P=332, THIS MAKES SENSE
            r = n0/n1
            if r>2 or r<0.5:
                print('stop here cause neutral is causing huge bending')


            direction, hasglitched = newdirection(poss[:,p], direction, r,p)

            if hasglitched ==1:
                glitchcount = glitchcount+1


            nextposs = poss[:,p] + ( interval * step * direction)#move along this arc for 1km  then recalc your position
            nextstraightposs = straightray[:,p] + + ( interval * step * constantdirection)
            straightray[:,p+1] = nextstraightposs
            poss[:,p+1] = nextposs #vsep doesnt allow strided arrays, so make a temp variable
            _,_, possalt[p+1] = spice.recgeo(nextposs, equatorialradii,flatteningcoefficient)
            _,_, straightpossalt[p+1] = spice.recgeo(nextstraightposs, equatorialradii,flatteningcoefficient)
        
            #need to have a catch for when the dtheta has been satisfied
            dtheta[p+1] = (spice.vsep( nextposs, MEX)) * (180 / math.pi)
            if dtheta[p] > angleseparation:
                break



        #shorten the arrays to just contain the non-zeros values
        dtheta = dtheta[dtheta != 0]
        
        possalt = possalt[possalt != 0]
        straightpossalt = straightpossalt[straightpossalt != 0]
        errors = np.zeros([1,p])
        straighterrors2 = np.zeros([1,p])
        straightrayalt = np.zeros([1,p])
        
        #ALT AND ANGLE PROGRESSION HAVE THE SAME LENGTH
        t=0
        alteration = alt[0]-possalt[0]
        for i in range(np.size(angleprogression,1)): 

            if angleprogression[0,i] > dtheta[t] and t<p :
                if t==0:
                    alteration = possalt[t]-alt[i]
                    straightalteration = straightpossalt[t] - alt[i]

                straightrayalt[0,t] = alt[i]
                
                bendalt = possalt[t]
                iteratedstraightrayalt = straightpossalt[t]
                error= straightrayalt[0,t] - bendalt +alteration
                straighterrors2[0,t] = straightrayalt[0,t] - iteratedstraightrayalt + straightalteration
                errors[0,t] = error
                t=t+1
            
            else:
                continue
            
        

        straighterrors = straight(ray,direction ,MEX, angleseparation,angleprogression, initialangle, xyzpoints,sv,interval)
        #if straighterrors2 == straighterrors, we can skip this

        straightcut = straighterrors[0, 0:p]
        #errors = errors[0]
        refractedcut = errors[0]
        # there is a glitch when the altitude becomes negative due to the spheriod not being considered. This is a patch to 
        # remove the point of inflexion. Further inspection shows that glitches occured in low altitude and not negative altidue, the expnential 
        # decay shape of the neutral refractivity profile led to small changes in alt, leading to huge changes in N, thus the refractive
        #  ratio can exceed two and the ray can bend extreamly, veering it of the straight ray and into the martian core, meaning it never touches the 2nd ionsphere
        
        #This shouldnt be required once a good starting angle is found, only required for graphing early ray propogations (rays with the largest misses)
        # for i in range(p-1): #
        #     if straightcut[i+1] < straightcut[i]:#inflexion
        #         straightcut[i+1] = straightcut[i] + (straightcut[i]-straightcut[i-1])
        #     if refractedcut[i+1] < refractedcut[i]:#inflexion
        #         refractedcut[i+1] = refractedcut[i] + (refractedcut[i]-refractedcut[i-1])


        results = straightcut - refractedcut
        dtheta = dtheta[0:-2]
        results = results[:-1]

        miss = results[-3] #take the last results, this is how many km away the ray is when dtheta is satisfied
        
        ploterrorvariation(dtheta[0:-2], results[0:-2], straightrayalt)

        iterationmax = iterationmax+1
   
def ploterrorvariation(dtheta, results, straightcut):
    
    smoothresults = UnivariateSpline(dtheta, results,s=0,k=4)
    secondderivative = smoothresults.derivative(n=2)

    #plt.plot(cut)
    #plt.plot(dtheta,possalt[0:-1] )
    #plt.plot(dtheta[0:-1],straightcut, label= 'straight' )
    #plt.plot(dtheta[0:-1],refractedcut, label = 'refracted' )

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    ax1.plot(dtheta, results)
    #ax1.legend()
    ax1.set_xlabel(r'$\delta \theta$') 
    ax1.set_ylabel('Refracted Propogation\'s Variation from Straight Line (km) ')
    ax2.set_ylabel('Variation\'\' ')
    plt.title(r'Variation between Refracted Path and Straight Line over $\delta \theta$')
    

    #ax2.plot(angleprogression[0,:],residual[0,:])
    ax2.plot(dtheta,secondderivative(dtheta))

    plt.show()

def newdirection(poss, direction, r,p):
    glitch = 0
    if r != 1: # if n0=n1 then dont change direction
        up =   poss / np.linalg.norm(poss) #vector going upwards radially if n0>n1
        ud = direction

        dotproduct = (up[0]*ud[0]) + (up[1]*ud[1]) + (up[2]*ud[2])
        if dotproduct<0: #if product is negative, reverse direction of unit possition (up)
            up = up * -1
            dotproduct = (up[0]*ud[0]) + (up[1]*ud[1]) + (up[2]*ud[2])
        subproduct = 1 - (r**2) * (1 - (dotproduct**2) )
        if subproduct <0:
            glitch = 1
            #do not alter direction
        else:
            c2 = math.sqrt(subproduct)
            direction = (r*ud) + (((r*dotproduct) - c2)*up)

    return direction, glitch



def straight( ray, initialdirection, MEX,  angleseparation, angleprogression, initialangle, xyzpoints,sv, interval):
    alt = ray[2,:]
   

    rawvector = xyzpoints[:,len(xyzpoints)-1]-xyzpoints[:,0]
    norm = np.linalg.norm(rawvector)
    rayunit = rawvector/norm
    interval =1


    marsrad = spice.bodvrd(sv.front, 'RADII', 3)
    flatteningcoefficient = ( marsrad[1][0] - marsrad[1][2] ) / marsrad[1][0]
    equatorialradii = marsrad[1][0]
    direction = initialdirection
    # everything must have a alt/ equiv dtheta
    poss =np.zeros([3,10000])
    dtheta = np.zeros([1,10000])
    possalt =np.zeros([1,10000])
    poss[:,0] = MEX
    step = 10  # this will scale down to 0.001 (1m)
    
    miniray = np.zeros([3,step])
    totalarcdistance = math.floor(alt[1] * angleseparation) # in km
    angleseparation = angleseparation * (180 / math.pi) #we want in degrees now
    #NOTE TO SELF, AT FIRST THE DIFFERNCE IN N MIGHT BE TO SMALL IN THE HIGH ALT. MAYBE JUMP TO WHERE THE IONO BEGINS THEN START THIS LOOP

    
    for p in range(100000//step): #10,000 is the max distance, this will never be reached
        if p == 0: #You need a before n and after n so for the begining go forwards two to get a n0 and n1
            point = poss[:,0] +0 # regeo needs non-strided arrays, this achieved by adding 0 (cheat)
            _,_, possalt[0,p] = spice.recgeo(point, equatorialradii,flatteningcoefficient)
            dtheta[0,p] = (spice.vsep( point, MEX))* (180 / math.pi) #should give 0 degrees
            for i in range(step):
                point = poss[:,0] + (interval * i * direction) #make a small 'step' long 3d vector
                miniray[:,i] = spice.recgeo(point, equatorialradii,flatteningcoefficient) #(lon; lat; alt)
            n0 = findrefractivity(miniray,step)
            
            
             
        else:
            n0=n1

        for i in range(step):
            point = poss[:,p] + (interval * i* direction) #make a small 'step' long 3d vector
            miniray[:,i] = spice.recgeo(point, equatorialradii,flatteningcoefficient)
        n1 = findrefractivity(miniray,step)
        

        #VECTOR CLACS [dont understand inuitivly, boiler plate] maybe this breaks if no r
        #SHOULD HIT THE IONO AT P=332, THIS MAKES SENSE
        r = n0/n1



        nextposs = poss[:,p] + (interval * step * direction)#move along this arc for 1km  then recalc your position
        poss[:,p+1] = nextposs #vsep doesnt allow strided arrays, so make a temp variable
        _,_, possalt[0,p+1] = spice.recgeo(nextposs, equatorialradii,flatteningcoefficient)
        #need to have a catch for when the dtheta has been satisfied
        dtheta[0,p+1] = (spice.vsep( nextposs, MEX))* (180 / math.pi)
        if dtheta[0,p] > angleseparation:
            break


    #shorten the arrays to just contain the non-zeros values
    dtheta = dtheta[dtheta != 0]
    
    possalt = possalt[possalt != 0]
    straighterrors = np.zeros([1,p])
    
    
    #ALT AND ANGLE PROGRESSION HAVE THE SAME LENGTH
    t=0
    alteration = alt[0]-possalt[0]
    for i in range(np.size(angleprogression,1)): #up tp 6478
        if t== 93:
            print('stophere')

        if angleprogression[0,i] > dtheta[t] and t<p :
            if t==0:
              alteration = possalt[t]-alt[i]

            straightrayalt = alt[i]
            bendalt = possalt[t]
            error= straightrayalt - bendalt +alteration
            straighterrors[0,t] = error
            t=t+1
        else:
            continue

    return straighterrors






#the issue must be here then, the function above makes sense.
def doppler(residual, totalperiods, dist): # the results of this is changing for the the number of iterations
    #dist needs to be altered, how many wavelenghts fit into the whole dist, that needs to be the iterable. wavelength is 61 cm. this is too small
    scale = 1 # scale =10, means we are itertating per 100 wavelenghts instead of 1000
    vacuumwavelength  = (constants.c / 437.1e6) /scale
    total1000periods = ((1/vacuumwavelength) * dist * scale) # to test if this is the rounding error in the geometricdistasnce 
    wavelength = np.ones([1,totalperiods])
    for i in range(totalperiods):
        wavelength[0,i] = (vacuumwavelength) * residual[0,i] #resisdual is only dist long

    electricdistance = np.sum(wavelength)
    geometricdistance = (vacuumwavelength*total1000periods)/scale # now this has not got a rounding error
    return electricdistance , geometricdistance
    
