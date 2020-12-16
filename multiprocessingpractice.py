import time
from multiprocessing import Pool, RawArray
import numpy as np
import csv
from os import path
import spiceypy as spice ; import spiceypy.utils.support_types as stypes
import pickle

import FindDopplerMain

def square(number):
    s=0
    for i in range (number):

        s += i*i
    return s



def mpsquare(numbers):
    starttime = time.time()
    p= Pool()
    results = p.map(square, numbers) #assign to the maximum number of cores on the machine 
    #print(results)
    p.close()
    p.join()
    endtime = time.time() - starttime 
    print(f"Processing {len(numbers)} numbers took {endtime} time using multiprocessing")

def nompsquare(numbers):
    starttime = time.time()
    results = []
    for i in numbers:
        results.append(square(i))
        endtime = time.time() - starttime 
    print(f"Processing {len(numbers)} numbers took {endtime} time using serial processing")

def ElectricDistance(numbers):
    
    starttime = time.time()
    p= Pool()
    #S = RawArray('d', 10)
    S = p.map(FindDopplerMain.ElectricCall,numbers) #assign to the maximum number of cores on the machine 
    print(S)
    p.close()
    p.join()
    endtime = time.time() - starttime 
    print(f"Processing {len(numbers)} numbers took {endtime} time using multiprocessing")
    return S

def GeometricDistance(numbers):
    
    starttime = time.time()
    p= Pool()
    #S = RawArray('d', 10)
    Geo = p.map(FindDopplerMain.GeoCall,numbers) #assign to the maximum number of cores on the machine 
    print(Geo)
    p.close()
    p.join()
    endtime = time.time() - starttime 
    print(f"Processing {len(numbers)} numbers took {endtime} time using multiprocessing")
    return Geo

#-----------------------------------------------------<VALUES TO EDIT REGULARLY>----------------------------------------
# If you only wish to analysis mutual [cross-link] occultation between MEX and TGO, then this is the only section that
# needs to be edited
start  = '2020 SEP 1'
stop   = '2020 SEP 3'
OCCSELECTION = 14 # Which occultation do you wish to see in Cosmographia? [optional]
here = path.abspath(path.dirname(__file__))
PathtoMetaKernel1 = here + '/TGO/mk/em16_plan.tm'
PathtoMetaKernel2 = here + '/MEX/mk/MEX_OPS.tm'
#PathtoMetaKernel3 = here + '/MarsRec_Odyssey/MRO_MO_mk.tm'
#-----------------------------------------------------------------------------------------------------------------------
#spice.furnsh(PathtoMetaKernel3)
spice.furnsh(PathtoMetaKernel1)
spice.furnsh(PathtoMetaKernel2)

        # epoch = spice.str2et(start)
        # print(epoch)


        # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~< calc the exact periods of egress occs
        # # Setting Variables
        # ingresslist = np.array([1.0] ,dtype = float)
        # egresslist = np.array([1.0] ,dtype = float)
        # etbeg = spice.str2et(start)
        # etend = spice.str2et(stop)

        # class SpiceVariables:
        #     obs  = '-41' # NAIF code for MEX
        #     target = '-143'# NAIF code for TGO ['EARTH'/'SUN'/ a groundstation etc]
        #     obsfrm = 'IAU_MARS'
        #     abcorr = 'NONE'
        #     crdsys = 'LATITUDINAL'
        #     coord  = 'LATITUDE'
        #     stepsz = 1.0 # Check every [300] seconds if there is an occultation
        #     MAXILV = 100000 #Max number of occultations that can be returned by gfoclt
        #     bshape = 'POINT'
        #     fshape = 'DSK/UNPRIORITIZED'
        #     front = 'MARS'
        #     fframe = 'IAU_MARS'
        #     TFMT = 'YYYY-MM-DD HR:MN:SC' # Format that Cosmographia understands

        # sv = SpiceVariables()



        # # Form a windows that gfoclt can populate
        # window = stypes.SPICEDOUBLE_CELL(2)
        # spice.wninsd(etbeg, etend,window)
        # occwindow = stypes.SPICEDOUBLE_CELL(sv.MAXILV)

        # #find occultation windows between the dates listed above [ most comp cost in this function]
        # spice.gfoclt('ANY',sv.front,sv.fshape,sv.fframe, sv.target, sv.bshape, '' , sv.abcorr, sv.obs, sv.stepsz, window, occwindow)

        # winsiz = spice.wncard( occwindow )# Find cardinality (number of windows)

        # #initialize lists to form dataframe
        # lon , lat, dist, sza, angle = ( np.ones(winsiz-1) for i in range(5))


        # # Enter the ingress epochs into a dataframe
        # occlist = np.ones((winsiz,3))
        # for i in range(winsiz):
        #     [ingress, egress] = spice.wnfetd(occwindow, i) # extract the begining and ends of the windows
        #     if i == 1 :
        #         ingresslist = ingress
        #         egresslist = egress
        #     else:
        #         ingresslist = np.append(ingresslist, [ingress])
        #         egresslist = np.append(egresslist, [egress])


        # exampleepoch = spice.timout(ingresslist[5], sv.TFMT)
        # print('et = ',ingresslist[5], 'Which in Julian is', exampleepoch)
        # with open('ingresslist.pkl', 'wb') as f:  # pickle the output due to multiprocessing not passing variables easily.
        #     pickle.dump([ingresslist], f)

#pickle the engress list 

if __name__ == '__main__':

    numbers = range(1000, 0, -1) # for now this number has to be kept under the logical core count of the computer (6), or you get repeated
                                # output, number of seconds you want to iterate prior to occ
    #numbers = range(5)
    #mpsquare(numbers)
    #nompsquare(numbers)
    S= ElectricDistance(numbers)
    Geo = GeometricDistance(numbers)
    Alt = FindDopplerMain.TangentPointAltitude(numbers)
    print('Combining electric distance and Real/Geometric Distance')
    S = np.vstack((S, Geo)) 
    S = np.vstack((S, Alt[:-1]))
    np.savetxt('ElectricDistance.csv', S, delimiter = ',')
    #np.savetxt('geo.csv', geo, delimiter = ',')
    print("this worked")

