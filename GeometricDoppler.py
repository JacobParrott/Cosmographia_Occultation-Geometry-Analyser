# a regular script for geometric doppler confinement will be required. include it in this project too
from os import path
import spiceypy as spice
from scipy import constants
import pickle
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

here = path.abspath(path.dirname(__file__))

PathtoMetaKernel1 = here + '/TGO/mk/em16_plan.tm'
PathtoMetaKernel2 = here + '/MEX/krns/mk/MEX_OPS.tm'

print(PathtoMetaKernel1)
print(PathtoMetaKernel2)

spice.furnsh(PathtoMetaKernel1)
spice.furnsh(PathtoMetaKernel2)



geometricdopplershift = np.zeros(600)
for time in range(599,0,-1):
    sc2scstates = spice.spkezr('-143', (659891734 - time), 'IAU_MARS', 'NONE', '-41')
    velocityvector = sc2scstates[0][3:6]
    velocityvector = velocityvector[0:3]
    positionalvector =  sc2scstates[0][0:3]
    positionalvector = positionalvector[0:3]
    velocityangle = spice.vsep( positionalvector, velocityvector) #rads
    relativevelocity = np.linalg.norm(velocityvector) * np.cos(velocityangle) 
    
    geometricdopplershift[time] = -(relativevelocity/constants.c) * 437.1e9 #conversion from km to m


#let's calculate doppler based on displacements and compare
displacement = np.zeros(600)
geometricdopplershift2 = np.zeros(600)
i=0
for time in range(599,0,-1):
    sc2scstates = spice.spkezr('-143', (659891734 - time), 'J2000', 'LT+S', '-41')
    positionalvector =  sc2scstates[0][0:3]
    positionalvector = positionalvector[0:3]
    displacement[i] = np.linalg.norm(positionalvector)
    velocity= displacement[i] - displacement[i-1] #change in displacment over 1 second
    geometricdopplershift2[time] =  - (velocity/constants.c) * 437.1e9

    i=i+1



plt.plot(geometricdopplershift[1:],'b',label = 'Built-in SPICE Relative Velocity Function')
plt.plot(geometricdopplershift2[1:-1],'r', label = 'Delta Displacement method')
plt.title('Both methods of finding the Geometric Doppler produce the same results')
plt.legend()
plt.show()


pd.DataFrame(geometricdopplershift).to_csv("geometricdopplershift.csv")
