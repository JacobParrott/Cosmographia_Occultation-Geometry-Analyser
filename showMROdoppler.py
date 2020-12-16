import csv ; import pandas as pd ; import numpy as np ; import matplotlib.pyplot as plt
from scipy import signal
import spiceypy as spice
import MRO_MO
import pickle
#.dat file has to be rewritten as csv file
# read flash.dat to a list of lists

#EXTRACTING THE ALT FROM THE FILE
MRO,MO, height, epoch = MRO_MO.extract_data()
alt = np.zeros(height-2)
for event in range(height-2):
    #calc the real distance per epoch
    mro = MRO[:,event] ; mo = MO[:,event]
    dis = mro-mo
    
    #calc the tangent point height
    sc2scvector = mro - mo
    displacement = np.linalg.norm(sc2scvector)
    sc2scunitvector = np.true_divide(sc2scvector, displacement)
    #marsrad = spice.bodvrd('MARS', 'RADII', 3)
    mro=mro+0
    _,alt[event] = spice.npedln(3396.19,3396.19, 3376.2, mro, sc2scunitvector)# this reports back 0 if it is below the ground




datContent = [i.strip().split() for i in open("./dopres_269.dat").readlines()] # this contains the resdoppler data

data = np.asarray(datContent)


#columns_to_keep = ['#time', 'b', 'c']
#df = pd.read_table("./dopres_253.dat", sep="\s+", usecols=columns_to_keep)
# write it as a new CSV file
        # with open("./flash.csv", "wb") as f:
        #     writer = csv.writer(f)
        #     writer.writerows(datContent)


        # def your_func(row):
        #     return row['x-momentum'] / row['mass']

        # columns_to_keep = ['#time', 'doppler']
        # dataframe = pd.read_csv("./flash.csv", usecols=columns_to_keep)
        # #dataframe['new_column'] = dataframe.apply(your_func, axis=1)

        # print(datafram)

doppler = np.zeros(np.size(data,axis=0))
time = np.zeros(np.size(data,axis=0))
for i in range(np.size(data,axis=0)):
    doppler[i] = float(data[i,1])

for i in range(np.size(data,axis=0)):
    time[i] = float(data[i,0])

#normalise time
#time = time - np.amin(time)
#epoch = epoch - np.amin(epoch)
#time = time + 242654465.18250176 #ET for midnight of that day
#datadelay = time[1] - epoch[1]
#print('residual doppler file begins at ', time[1], '\nAnd Geometry file begins at', epoch[1], '\nDoppler file is is ',time[1]-epoch[1], 'seconds later')
#time = time + datadelay

#doppler = signal.savgol_filter(doppler,101,3)

#applying a butterworth filter to avoid shot noise in rotation error
# period = len(time)
# maxfreq = 100
# nyquist = 0.5
# num,denom = signal.butter(2, 0.02, btype = 'low', analog = False)
# doppler = signal.filtfilt(num, denom, doppler)

with open('MRO_MO_doppler269.pkl', 'wb') as f:  # pickle the output due to multiprocessing not passing variables easily.
    pickle.dump([doppler, time], f)

fig, ax = plt.subplots(1,1)
ax.plot(time,doppler)
ax.set_title('Residual Doppler from MRO -> Odssey on 10/09/2007')
ax.set_xlabel('Time(s)[end == occultation epoch]')
ax.set_ylabel('Residual Doppler (Hz)')
#ax2 = ax.twinx()
#ax2.plot(epoch[ 0:-4], alt)
plt.show()