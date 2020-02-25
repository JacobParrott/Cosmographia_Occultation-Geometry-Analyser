import numpy as np
import spiceypy as spice
import math
import JSONtemplate as JS
import pandas as pd
import json
from PIL import Image
import matplotlib.pyplot as plt
import cartopy.crs as ccrs #import the coordinate refernece system


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
        #template['items'][itemcounter]['geometry']['emitters'][0]['startTime'] = spice.timout((et - i), sv.TFMT)
        #template['items'][itemcounter]['geometry']['emitters'][0]['endTime'] = endtimeFMT
        template['items'][itemcounter]['geometry']['emitters'][0]['velocityVariation'] = (profilespread + 1)
        template['items'][itemcounter]['geometry']['emitters'][0]['endSize'] = profileintensity
        itemcounter = itemcounter + 1

    # serialise the formed profile into a .json that cosmographia can read
    with open(' Profile.json', 'w') as file:
        json.dump(template, file, indent=3)


# Find the location of the lowest point of the occultation
def Location(et, sv):
    Coords = np.ones(3)
    [tgopos, _] = spice.spkpos(sv.front, et, sv.fframe, 'NONE', sv.target)
    [sc2scvector,_] = spice.spkpos(sv.target, et, sv.fframe, 'NONE', sv.obs)
    displacement = np.linalg.norm(sc2scvector)
    sc2scunitvector = np.true_divide(sc2scvector, displacement)
    marsrad = spice.bodvrd(sv.front, 'RADII', 3)  # Extract the triaxial dimensions of Mars
    # For the ray that connects MEX and TGO, find the point on this ray that is closest to the Martian surface
    [nearestpoint,_] = spice.npedln(marsrad[1][0], marsrad[1][1], marsrad[1][2], tgopos, sc2scunitvector)
    [_, lon, lat] = spice.reclat(nearestpoint)
    lon = lon * (-180 / math.pi) # Rad -> Deg , frame inversion required
    lat = lat * (-180 / math.pi)
    #Coords[0] = lon
    #Coords[1] = lat
    return lon,lat, displacement, nearestpoint



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
    path_to_pic = file_location + '/Pictures/2k_mars.jpg'
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