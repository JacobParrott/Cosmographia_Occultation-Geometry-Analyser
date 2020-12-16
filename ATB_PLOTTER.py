#PLOT THE ATB DATA WITH A TRUNCATION AND 10 HZ SMOOTHING
from scipy import signal ; from scipy import interpolate ;import pandas as pd ; import numpy as np ; import spiceypy as spice ;  import pickle ; import matplotlib.pyplot as plt

# order of columns :# 0]Date Time 2] Phase(rad)   3] Loopstate 4]Lvl(dB) 5] EsNo    6] Timeofs(s) 6] Lock 6]XpdrAGC   7] Freq(Hz)    8]C/No

df = pd.read_fwf('muttocc_C.txt', sep =  ' ', header =None )
df.columns = ['date', 'time' , 'phase'  ,'loopstate' , 'level' , 'spectralefficiency' , 'timeoffset' , 'lock' , 'AGC' , 'freq' ,'CarrierNoise']
#df.columns = ['EphemerisTime','mexlat','mexlon','mexalt','tgolat','tgolon','tgoalt','DShift']
print(df.head())

freq = df['freq'].to_numpy()
geodoppler= np.genfromtxt('geometricdopplershift.csv', delimiter=',')

geodoppler = geodoppler[2:,1] #599 units
# convert column 'freq' to an array

level = df['level'].to_numpy()
Noiseratio = df['CarrierNoise'].to_numpy()

#mylist = df.index[df['loopstate'] == 'Lock'].tolist()

# filtered = df[df.loopstate.isin(['Locked'])]
# level = filtered['level'].to_numpy()
# Noiseratio = filtered['CarrierNoise'].to_numpy()


t=0

lockstate = df['loopstate'].tolist() #using the loopstate to assist with filtering
for i in range(1,(len(lockstate)-1)):
    if lockstate[i+1] == 'Unlocked':
        freq[i] = freq[t]
    else:
        if lockstate[i-3] == 'Locked':
            t=i

xnew = np.arange(0,600,0.1)
freq = signal.resample(freq, 6000) #10 hz over 10 mins
#geodoppler = signal.resample(geodoppler, 6000) #10 hz over 10 mins
freq = signal.medfilt(freq ,61)
freq = signal.savgol_filter(freq ,1001,3)
freq = signal.savgol_filter(freq ,501,3)




#make a for loop for filtering 

t=0

for i in range(1,(len(lockstate)-1)):
    if lockstate[i+1] == 'Unlocked':
        Noiseratio[i] = Noiseratio[t]
    else:
        if lockstate[i-3] == 'Locked':
            t=i
xnew = np.arange(0,600,0.1)
Noiseratio = signal.resample(Noiseratio, 6000) #10 hz over 10 mins
Noiseratio = signal.medfilt(Noiseratio ,61)
Noiseratio = signal.savgol_filter(Noiseratio ,1001,3)
Noiseratio = signal.savgol_filter(Noiseratio ,501,3)

                    # freq = signal.savgol_filter(freq ,101,3)
                    # #filter visually if the signal exeeds bounds of 1%
                    # for i in range(60,(len(lockstate)-1)): #move from front to back
                    #     avgfreq = np.mean(freq[i-3:i])
                    #     if (freq[i]> 1.01*avgfreq) or (freq[i]< 0.99*avgfreq):
                    #         freq[i] = avgfreq



#need to filter out the spikes, find the index of the 'unlock' states
#use those indexs to set the freq values with the same index to a previos value before the 'unlock'
#this way you dont remove any time, which will lead to a drift over time
                # i = 0
                # window_size = 20
                # moving_averages = []
                # while i < len(resampled) - window_size + 1:
                #         this_window = resampled[i : i + window_size]
                #         window_average = sum(this_window) / window_size
                #         moving_averages.append(window_average)
                #         i += 1
                # resampled = moving_averages

                                    # 

#need to inspect the filtering process of the new dataframe
#filtered.to_csv(r'filtered.csv')

        # num,denom = signal.butter(4,0.5, btype = 'low', analog = False)
        # resampled = signal.filtfilt(num, denom, resampled)




                        # resampled = signal.savgol_filter(resampled ,101,3)


                # i = 0
                # window_size = 20
                # moving_averages = []
                # while i < len(resampled) - window_size + 1:
                #         this_window = resampled[i : i + window_size]
                #         window_average = sum(this_window) / window_size
                #         moving_averages.append(window_average)
                #         i += 1
                # resampled = moving_averages



                                    # resampled = signal.savgol_filter(resampled ,401,3)

# 

# f = interpolate.interp1d(xnew,resampled ,kind = 'cubic')
# #f.set_smoothing_factor(10)
# xnew2 = np.arange(5, 600,0.01)
# interped = f(xnew2)



# dfstates = pd.read_fwf('TableOfStates.txt', sep =  ' ', header =None )
# dfstates.columns = ['EphemerisTime','mexlat','mexlon','mexalt','tgolat','tgolon','tgoalt','DShift']
# dfstates['atmos'] = interped
#np.savetxt('freq.csv', interped, delimiter = ',')

fig, ax = plt.subplots(1)
plt.plot(xnew[0:5000], freq[0:5000], 'b', label ='Measured Data') 
plt.title('Measured Carrier Offset and SPICE simulated Geometric Doppler [truncated at Occultation Epoch]')
plt.ylabel('Carrier Offset [Hz]')
#ax1 = ax[0].twinx()

plt.plot(geodoppler[0:500], 'r', label = 'SPICE Simulated Data') 
#ax1.set_title('SPICE Simulated Geometric Doppler Shift')
# ax1.set_ylabel('Doppler Shift [Hz]')
# #plt.plot(xnew,resampled) #[0:-(window_size-1)]
# #plt.plot(xnew2,interped, 'r')
plt.xlabel('time [s]')
plt.legend()


plt.show()
