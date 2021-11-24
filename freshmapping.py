import pandas as pd
from os import path
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
from cartopy.feature.nightshade import Nightshade
from datetime import datetime
import numpy as np
import shapely.geometry as sgeom
from scipy.ndimage.filters import gaussian_filter
from math import radians, atan2, sin, cos, sqrt
from matplotlib.image import imread
import spiceypy as spice

import main


def GlobePlotter(solarlat, solarlon, lattrace, lontrace):
    midlat = np.mean(lattrace)
    midlon = np.mean(lontrace)
    projection = ccrs.Orthographic(midlon, midlat)
    ax = plt.axes(projection=projection)
    file_location = path.abspath(path.dirname(__file__))
    path_to_pic = file_location + '/images/2k_mars.jpg'
    source_proj = ccrs.PlateCarree()
    ax.imshow(imread(path_to_pic), origin='upper', transform=source_proj,
              extent=[-180, 180, -90, 90])

    # plot the ground feature labels
    df = pd.read_csv('mars.csv', encoding='latin-1')
    df.sort_values(by=["Diameter"], inplace=True, ascending=False)
    minlen = 30    # plot only the largest 30 objects on the surface of mars. Maybe update to the the best known
    df = df.head(min(minlen, len(df)))
    for index, row in df.iterrows():
        text = row['Clean_Feature_Name']
        x, y, s = row['Center_Longitude'], row['Center_Latitude'], row['Diameter']

        ax.text(x, y, text, transform=ccrs.PlateCarree(),
                ha='left', va='center', fontsize=8, color='#ebc334')
        ax.scatter(x, y, transform=ccrs.PlateCarree(),
                   s=10, color='#ebc334', edgecolor='None', lw=0)

    # adjusting the gridlines to fit the occtrace location [10 MIGHT BE TO COARSE]
    # STILL NEED BETTER SOLUTION
    minlat = (np.floor(np.min(lattrace)/10))*10
    maxlat = (np.ceil(np.max(lattrace)/10))*10
    minlon = (np.floor(np.min(lontrace)/10))*10
    maxlon = (np.ceil(np.max(lontrace)/10))*10
    track = sgeom.LineString(zip(lontrace, lattrace))
    ax.add_geometries([track], source_proj,
                      edgecolor='#C852C8', linewidth=8, alpha=0.5, facecolor='none')
    gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=2,
                      color='k', alpha=0.2, linestyle='--', draw_labels=True)
    gl.top_labels = False
    gl.left_labels = False
    gl.right_labels = False
    gl.xlines = True
    gl.ylocator = mticker.FixedLocator(
        list(np.arange(int(minlat), int(maxlat), 10)))
    gl.xlocator = mticker.FixedLocator(
        list(np.arange(int(minlon), int(maxlon), 10)))
    gl.xlabel_style = {'color': 'k'}

    # <<<< this involves a custom nightshade function, go into the cartopy libs and edit to include lat lon
    ax.add_feature(Nightshade(solarlat, solarlon))

    # plot occultation track
    ax.add_geometries([track], source_proj,
                      edgecolor='#C852C8', linewidth=8, alpha=0.5, facecolor='none')
    plt.show()


def f():

    here = path.abspath(path.dirname(__file__))

    PathtoMetaKernel1 = here[:-39] + \
        '/git/exomars2016/kernels/mk/em16_ops.tm'
    PathtoMetaKernel2 = here[:-39] + \
        '/git/mars-express/kernels/mk/MEX_OPS.tm'
    print(PathtoMetaKernel1)
    print(PathtoMetaKernel2)
    spice.furnsh(PathtoMetaKernel1)
    spice.furnsh(PathtoMetaKernel2)

    sv = main.SpiceVariables()
    et = 657605100
    # get an actual occultation trace
    [trace, altTrace] = main.occSurfaceTrace(et, sv)
    #lattrace = np.linspace(10, 40, 10)
    #lontrace = np.linspace(90, 170, 10)
#    print('the starting locations are ', lontrace[1] + ',', lattrace[1])

    lattrace = trace[:, 0]
    lontrace = trace[:, 1]
    print('the starting locations are ', lontrace[1], ',', lattrace[1])

    GlobePlotter(20, 10, lattrace, lontrace)


if __name__ == '__main__':
    f()
