#THIS SCRIPT IS ONYL REQUIRED IF YOU CANNOT IMPORT THE ACTUAL SPICE KERNALS FOR MRO AND MO

import pandas as pd ; import numpy as np ; import spiceypy as spice ; from os import path ; from multiprocessing import Pool, RawArray ; import pickle
import FindDopplerMain
import main


def MRO_MO(event):
    with open('MRO_MO_states.pkl', 'rb') as f:  # Python 3: open(..., 'rb')
        MRO, MO, height, et = pickle.load(f)
    print('progress = ', format( ((event/height)*100), '.3f') )
    mro = MRO[:,event] ; mo = MO[:,event] ; time = et[event]
    referencedirection = [0,0,0]
    dis = mro-mo
    Distance = np.linalg.norm(dis)
    initialangle,xyzpoints= FindDopplerMain.producegeometrymeter(mro,mo, time) #make 2D for huge speed improvements.
    Bending, S , referencedirection = FindDopplerMain.flatbending(xyzpoints,initialangle, mro,mo, referencedirection)
    return S

def extract_data():
    text_file = open('ody_to_mro_2007-253.txt', 'r')
    lines = text_file.read().split('\n')
    print(lines[10])
    height = np.size(lines, axis=0)
    states = np.zeros((height, 17))
    year , doy , hr , minute , sec , x1 , y1 , z1 , x2 , y2 , z2 , vx1 , vy1 , vz1 , vx2 , vy2 , vz2, epoch = (([0]*(height+2)) for i in range(18))
    for i in range(height-1):
        state = lines[i].split(' ')
        year[i] = state[0]
        doy[i] = state[1]
        hr[i] = state[2]
        minute[i] = state[3]
        sec[i] = state[4]
        x1[i] = float(state[5]) ; y1[i] = float(state[6]) ; z1[i] = float(state[7])
        x2[i] = float(state[8]) ; y2[i] = float(state[9]) ; z2[i] = float(state[10])

        vx1[i] = float(state[11]) ; vy1[i] = float(state[12]) ; vz1[i] = float(state[13])
        vx2[i] = float(state[14]) ; vy2[i] = float(state[15]) ; vz2[i] = float(state[16])

        stringepoch = year[i] + '-' +  doy[i] + '::' + hr[i] + ':' + minute[i] + ':' + sec[i] #convert times into SPICE comprehendable string
        epoch[i] = spice.str2et(stringepoch)
        #'YYYY DOY HR MN SC'

    MRO = np.array([x1,y1,z1])
    MO = np.array([x2,y2,z2])
    et = np.array(epoch)
    
    with open('MRO_MO_states.pkl', 'wb') as f:  # pickle the output due to multiprocessing not passing variables easily.
        pickle.dump([MRO, MO, height, et], f)

    return MRO, MO, height, epoch


#reading a large ASCII [separated by spaces]
##year doy hr min sec x1[km] y1[km] z1[km] x2[km] y2[km] z2[km] vx1[km/s] vy1[km/s] vz1[km/s] vx2[km/s] vy2[km/s] vz2[km/s] 





#IMPORT PREV KERNELS TO GET LEAP SECONDS AND MARS RADII
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
here = path.abspath(path.dirname(__file__))

PathtoMetaKernel1 = here + '/TGO/mk/em16_plan.tm'
PathtoMetaKernel2 = here + '/MEX/krns/mk/MEX_OPS.tm'
PathtoMetaKernel3 = here + '/MarsRec_Odyssey/MRO_MO_mk.tm'

spice.furnsh(PathtoMetaKernel1)
spice.furnsh(PathtoMetaKernel2)
spice.furnsh(PathtoMetaKernel3)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




if __name__ == '__main__':




    MRO,MO, height,epoch = extract_data()
    #height = 12 # just to test the core and speed things up
    geo = np.zeros(height-2)
    alt = np.zeros(height-2)
    
    for event in range(height-2):
        #calc the real distance per epoch
        mro = MRO[:,event] ; mo = MO[:,event]
        dis = mro-mo
        geo[event] = np.linalg.norm(dis)
        #calc the tangent point height
        sc2scvector = mro - mo
        displacement = np.linalg.norm(sc2scvector)
        sc2scunitvector = np.true_divide(sc2scvector, displacement)
        #marsrad = spice.bodvrd('MARS', 'RADII', 3)
        mro=mro+0
        _,alt[event] = spice.npedln(3396.19,3396.19, 3376.2, mro, sc2scunitvector)# this reports back 0 if it is below the ground


    #test in a single instance, !!then implement multitasking!!
    mro = MRO[:,1] ; mo = MO[:,1]


    
    events = range(height-2)
    p = Pool()
    S = p.map(MRO_MO ,events)
    print(S)
    p.close()
    p.join()

    S = np.vstack((S, geo))
    S = np.vstack((S, alt))


    np.savetxt('MRO_MO_ElectricDistance.csv', S, delimiter = ',')

    print('Electric Distance = ',S)
    print( 'Real Distance = ', geo)
    print('Tangent point alt = ', alt)
    print('Execution Successful')


