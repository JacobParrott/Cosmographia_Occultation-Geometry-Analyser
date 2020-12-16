import numpy as np
import matplotlib.pyplot as plt 
data = np.genfromtxt("DataATB.txt", dtype =None)
print(data[1])

print(data[1][9])

Freqs = []
i=0

Freqs = []
i=0
while i < 6031:
    Freq = data[i][9]
    Freqs = np.append(Freqs, Freq)
    
    i=i+1
print(Freqs)
len(data)


plt.style.use('dark_background')
plt.plot(Freqs)
plt.show()