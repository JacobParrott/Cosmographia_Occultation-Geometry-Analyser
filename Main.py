
import numpy as np
import spiceypy as spice
import spiceypy.utils.support_types as stypes
import math
import pandas as pd
import json

#Import custom modules
import profileformer
import JSONtemplate as JS


#-----------------------------------------------------<VALUES TO EDIT REGULARLY>----------------------------------------
# If you only wish to analysis mutual [cross-link] occultation between MEX and TGO, then this is the only section that
# needs to be edited
start  = '2020 JAN 1'
stop   = '2020 JAN 6'
OCCSELECTION = 18 # Which occultation do you wish to see?
PathtoMetaKernel1 = '/home/jparrott/scratch/jparrott/windows_mocc/tgo_kerns/mk/em16_plan.tm'
PathtoMetaKernel2 = '/home/jparrott/scratch/jparrott/windows_mocc/mex_kerns/mk/MEX_OPS.tm' #[ ommit if only one mk]
#-----------------------------------------------------------------------------------------------------------------------


spice.furnsh(PathtoMetaKernel1)
spice.furnsh(PathtoMetaKernel2)
obs  = '-41' # NAIF code for MEX
target = '-143'# NAIF code for TGO ['EARTH'/'SUN'/ a groundstation etc]
obsfrm = 'IAU_MARS'
abcorr = 'NONE'
crdsys = 'LATITUDINAL'
coord  = 'LATITUDE'
stepsz = 300.0 # Check every 300 seconds if there is an occultation
MAXILV = 100000 #Max number of occultations that can be returned by gfoclt
bshape = 'POINT'
fshape = 'DSK/UNPRIORITIZED'
front = 'MARS'
fframe = 'IAU_MARS'
TFMT = 'YYYY-MM-DD HR:MN:SC' # Format that Cosmographia understands


# Setting Variables
ingresslist = np.array([1.0] ,dtype = float)
etbeg = spice.str2et( start )
etend = spice.str2et( stop  )

# Form a windows that gfoclt can populate
window = stypes.SPICEDOUBLE_CELL(2)
spice.wninsd(etbeg, etend,window)
occwindow = stypes.SPICEDOUBLE_CELL(MAXILV)

#find occultation windows between the dates listed above [ most comp cost in this function]
spice.gfoclt('ANY',front,fshape,fframe, target, bshape, 'J2000' , abcorr, obs, stepsz, window, occwindow)

winsiz = spice.wncard( occwindow )# Find cardinality (number of windows)

# Inter the ingress epochs into a dataframe
occlist = np.ones((winsiz,3))
for i in range(winsiz):
    [ingress, eegress] = spice.wnfetd(occwindow, i) # extract the begining and ends of the windows
    if i == 1 :
        ingresslist = ingress
    else:
        ingresslist = np.append(ingresslist, [ingress])
occ=np.transpose(range(winsiz-1))
df = pd.DataFrame(ingresslist, columns=['Time'])

# Select and occultation of interest, calculate the shape of it's profile, then add this to a dateframe
et = df.Time[OCCSELECTION]
result = profileformer.occgradient(front,et , fframe,obs,target )
profile = result[0]
lowesttime = result[1]
highesttime = result[2]
lowesttimeFMT = spice.timout((et-result[1]),TFMT)
highesttimeFMT = spice.timout((et-result[2]),TFMT)
endtimeFMT = spice.timout(et,TFMT)
profile = np.transpose(profile)
profiledataframe = pd.DataFrame(profile[:][0:highesttime], columns = ['X', 'Y', 'Z'])

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
vs = vs.repeat(200,axis=0)
vt = vs[:][0:highesttime]
velocitiesdataframe = pd.DataFrame(vt, columns = ['dX', 'dY', 'dZ'])
finalprofile = pd.concat([profiledataframe,velocitiesdataframe],axis=1)


# Construct a JSON template depending on the size of the profile, split into sections of 10 points [smaller sections,
# more smoothly the profile forms over time]
blockcount = math.floor(highesttime/10)*10
i=0
while i < (blockcount/10):
  JS.JSONiterated = JS.JSONiterated +JS.JSONiterable
  i=i+1
JSONtemplate = JS.JSONstart + JS.JSONiterated + JS.JSONend


template = json.loads(JSONtemplate) #load the template created above  so it can be editted as JSON

itemcounter = 0

#convert states into single string for cosmographia to understand [cosmogrpahia computes as
# 3 postions and 3 velocities, no need for '\n']
for i in range(0,blockcount,10):
  block = (finalprofile[i:i+10].to_numpy())*(-1) #inverse due to reference frame inversion
  ProfileBlock = block.tolist() # convert to list for the JSON format

  states=[]
  for p in range(10):
    states.extend(ProfileBlock[p]) #turn states into one long line to remove the '\n', which confuses cosmographia

    # vary the spread and intensity of the profile as it increases with alitude and the atmosphere becomes less dense
  profilespread = math.floor(itemcounter * (5/(blockcount/10)))
  profileintensity = math.floor(itemcounter * (30 / (blockcount / 10)))

  # edit the large JSON file iterativly, with 'itemcounter'moving down the JSON items (blocks of 10 profile points).
  # adding the states(position and velocity of the polygons), name of item (so comsmographia doesnt overwrite them),
  # the start time, end time(epoch of occultation), the spread and finally the intensity of colour
  template['items'][itemcounter]['geometry']['emitters'][0]['generator']['states'] = states
  template['items'][itemcounter]['name'] = 'Profile ' + str(itemcounter)
  template['items'][itemcounter]['geometry']['emitters'][0]['startTime']  = spice.timout((et-i),TFMT)
  template['items'][itemcounter]['geometry']['emitters'][0]['endTime'] = endtimeFMT
  template['items'][itemcounter]['geometry']['emitters'][0]['velocityVariation'] = (profilespread +1)
  template['items'][itemcounter]['geometry']['emitters'][0]['endSize'] = profileintensity
  itemcounter=itemcounter+1

#serialise the formed profile into a json that cosmographia can read
with open( ' Profile.json'  , 'w') as file:
    json.dump(template,file, indent = 3)










