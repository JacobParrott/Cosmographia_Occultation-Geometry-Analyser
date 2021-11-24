# mapping invoker
import numpy as np
import spiceypy as spice
import spiceypy.utils.support_types as stypes
import pandas as pd
from os import path
import matplotlib.pyplot as plt
import time
from tqdm import tqdm
import csv
from multiprocessing import Pool
import math
from scipy import constants
from PIL import Image
import cartopy.crs as ccrs  # import the coordinate refernece system


# Find the location of the lowest point of the occultation


def Location(et, ingress, sv, when):
    Coords = np.ones(3)
    [tgopos, _] = spice.spkpos(sv.front, et-when, sv.fframe, 'NONE', sv.target)
    [mexpos, _] = spice.spkpos(sv.front, et-when, sv.fframe, 'NONE', sv.obs)
    [states, _] = spice.spkezr(sv.target, et-when, sv.fframe, 'NONE', sv.obs)
    sc2scvector = states[0:3]
    velocity = states[3:6]
    relativespeed = np.linalg.norm(velocity)
    # e9 because we are converting from km to m (SPICE outputs km, but constants in m)
    veldopp = (relativespeed/constants.c) * 437.1e9
    displacement = np.linalg.norm(sc2scvector)
    sc2scunitvector = np.true_divide(sc2scvector, displacement)
    # Extract the triaxial dimensions of Mars
    marsrad = spice.bodvrd(sv.front, 'RADII', 3)
    # For the ray that connects MEX and TGO, find the point on this ray that is closest to the Martian surface
    [nearestpoint, alt] = spice.npedln(
        marsrad[1][0], marsrad[1][1], marsrad[1][2], tgopos, sc2scunitvector)
    # THERE IS MORE SETTINGS ON THIS
    [radius, lon, lat] = spice.reclat(nearestpoint)
    # Rad -> Deg , frame inversion required (hence the negative 180)
    lon = 180 - (lon * (-180 / math.pi))
    lat = lat * (-180 / math.pi)

    MexNadirTGOAngle = spice.vsep(-mexpos, -sc2scvector)
    MexNadirTGOAngle = MexNadirTGOAngle * (180/math.pi)

    # produce a string of the date and time, because an ephemeris time is not human-readable
    date_time = spice.timout(et, 'MM-DD HR:MN:SC')
    ingress_date_time = spice.timout(ingress, 'MM-DD HR:MN:SC')
    return lon, lat, displacement, nearestpoint, alt, relativespeed, date_time, ingress_date_time, veldopp, MexNadirTGOAngle


def charter(lon, lat, beg, stop, file_location):
    path_to_pic = file_location + '/images/2k_mars.jpg'
    raw_image = Image.open(path_to_pic)
    img = np.asarray(raw_image)  # convert to array
    globe = ccrs.Globe(semimajor_axis=285000.,
                      semiminor_axis=229000., ellipse=None)
    crs = ccrs.PlateCarree(globe=globe)
    extent = (-895353.906273091, 895353.906273091, 447676.9531365455, -
              447676.9531365455)  # adjustments of image for Robinson projection
    projection = ccrs.Robinson()
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=projection)
    ax.imshow(raw_image, transform=crs, extent=extent)

    title = 'Occultation locations between '+beg[5:]+' and ' + stop[5:]
    plt.title(title)
    ax.plot(lon, lat, 'o', c='#bef9b9', transform=ccrs.PlateCarree())
    plt.show()

# EXTENT MUST HAVE THE SAME DIMENTIONS AS GLOBE


def Newcharter(lon, lat, beg, stop, file_location):
    file_location = path.abspath(path.dirname(__file__))
    path_to_pic = file_location + '/images/2k_mars.jpg'
    raw_image = Image.open(path_to_pic)
    img = np.asarray(raw_image)  # convert to array
    globe = ccrs.Globe(semimajor_axis=285,
                       semiminor_axis=285, ellipse=None)
    reference = ccrs.Mollweide()
    extent = (-895353.906273091, 895353.906273091, 447676.9531365455, -
              447676.9531365455)  # adjustments of image for Robinson projection
    projection = reference

    fig = plt.figure()
    ax = plt.axes(projection=projection)
    ax.imshow(raw_image, transform=reference, extent=extent)

    plt.show()

# 🎇from https://scitools.org.uk/cartopy/docs/v0.15/examples/aurora_forecast.html

# def fill_dark_side(ax, time=None, *args, **kwargs):
#     """
#     Plot a fill on the dark side of the planet (without refraction).

#     Parameters
#     ----------
#         ax : matplotlib axes
#             The axes to plot on.
#         time : datetime
#             The time to calculate terminator for. Defaults to datetime.utcnow()
#         **kwargs :
#             Passed on to Matplotlib's ax.fill()

#     """
#     lat, lng = sun_pos(time)
#     pole_lng = lng
#     if lat > 0:
#         pole_lat = -90 + lat
#         central_rot_lng = 180
#     else:
#         pole_lat = 90 + lat
#         central_rot_lng = 0

#     rotated_pole = ccrs.RotatedPole(pole_latitude=pole_lat,
#                                     pole_longitude=pole_lng,
#                                     central_rotated_longitude=central_rot_lng)

#     x = np.empty(360)
#     y = np.empty(360)
#     x[:180] = -90
#     y[:180] = np.arange(-90, 90.)
#     x[180:] = 90
#     y[180:] = np.arange(90, -90., -1)

#     ax.fill(x, y, transform=rotated_pole, **kwargs)


class SpiceVariables:
    obs = '-41'  # NAIF code for MEX (-41)
    # NAIF code for TGO (-143)['EARTH'/'SUN'/ a groundstation etc]
    target = '-143'
    obsfrm = 'IAU_MARS'
    abcorr = 'NONE'
    crdsys = 'LATITUDINAL'
    coord = 'LATITUDE'
    stepsz = 2.0  # Check every 2 seconds if there is an occultation
    MAXILV = 100000  # Max number of occultations that can be returned by gfoclt
    bshape = 'POINT'  # Rx shape
    fshape = 'ELLIPSOID'
    front = 'MARS'
    fframe = 'IAU_MARS'
    TFMT = 'YYYY-MM-DD HR:MN:SC'  # Format that Cosmographia understands


start = '2020 NOV 30'
stop = '2020 DEC 2'
# Which occultation do you wish to see in Cosmographia? [optional]
OCCSELECTION = 2
here = path.abspath(path.dirname(__file__))

PathtoMetaKernel1 = 'C:/Users/Jacob/Documents/Doppler-Simulation-for-Mutual-Occultation/TGO/mk/em16_ops.tm'
PathtoMetaKernel2 = 'C:/Users/Jacob/Documents/Doppler-Simulation-for-Mutual-Occultation/MEX/mk/MEX_OPS.tm'

print(PathtoMetaKernel1)
print(PathtoMetaKernel2)

spice.furnsh(PathtoMetaKernel1)
spice.furnsh(PathtoMetaKernel2)
# spice.furnsh(PathtoMetaKernel3)

sv = SpiceVariables()


# Setting Variables
ingresslist = np.array([1.0], dtype=float)
egresslist = np.array([1.0], dtype=float)
etbeg = spice.str2et(start)
etend = spice.str2et(stop)


# Form a windows that gfoclt can populate
window = stypes.SPICEDOUBLE_CELL(2)
spice.wninsd(etbeg, etend, window)
occwindow = stypes.SPICEDOUBLE_CELL(sv.MAXILV)

# find occultation windows between the dates listed above [ most comp cost in this function]
spice.gfoclt('ANY', sv.front, sv.fshape, sv.fframe, sv.target,
             sv.bshape, '', sv.abcorr, sv.obs, sv.stepsz, window, occwindow)

winsiz = spice.wncard(occwindow)  # Find cardinality (number of windows)

# initialize lists to form dataframe
lon, lat, dist, sza, angle = (np.ones(winsiz-1) for i in range(5))


# Enter the ingress epochs into a dataframe
occlist = np.ones((winsiz, 3))
for i in range(winsiz):
    # extract the begining and ends of the windows
    [ingress, egress] = spice.wnfetd(occwindow, i)
    if i == 1:
        ingresslist = ingress
        egresslist = egress
    else:
        ingresslist = np.append(ingresslist, [ingress])
        egresslist = np.append(egresslist, [egress])

# form the dataframe
occs = pd.DataFrame(ingresslist, columns=['Time'])
occs['Ingress'] = ingresslist
occs['Egress'] = egresslist

date_time = [""]*(winsiz-1)
ingress_date_time = [""]*(winsiz-1)
localtime = [""]*(winsiz-1)
lon = np.ones([winsiz-1, 1])
veldopp = np.ones([winsiz-1, 1])
lat = np.ones([winsiz-1, 1])
dist = np.ones([winsiz-1, 1])
sza = np.ones([winsiz-1, 1])
clearanceangle = np.ones([winsiz-1, 1])
Rx = np.ones([winsiz-1, 1])
MexNadirTGOAngle = np.ones([winsiz-1, 1])
# grazingangle = np.ones([winsiz-1,1]) # HERE

for i in tqdm(range(winsiz-1)):

    # try to plot the location on a map with cartopy
    lon[i], lat[i], dist[i], nearestpoint, alt, speed, date_time[i], ingress_date_time[i], veldopp[i], MexNadirTGOAngle[i] = Location(
        occs.Time[i], occs.Ingress[i], sv, 0)  # FUNCTION NEEDS NEW NAME


# plot all of the tangent points onto the surface of mars
Newcharter(lon, lat, start, stop, here)
