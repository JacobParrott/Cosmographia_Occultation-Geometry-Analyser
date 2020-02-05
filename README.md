# Occultation Geometry Analyser
![](https://github.com/JacobParrott/OccultationProfiler/blob/master/images/Coverimage.png)
## What this project aims to achieve
In Short - It aims to create a tool to analyse occultation geometry to aid in observation planning. This tool will be built using the Spiceypy package so that the SPICE toolset( developed by NAIF ) can be employed. The first few releases of this project will aim to produce an automated Cosmographia catalogue to visualise the shape of an occultation profile ( as shown above ). Later on this project will be updated to produce further details about individual occultations, such as:
+ Precise coordinates ( with map )
+ Solar Zenith Angles
+ Local times
+ Distance between target and observer
+ Angle from lowest tangent point to highest point
+ Expected Rx (dB) 

>Finally, with all the above calculated, occultations with the **perfect characteristics** can be selected.

*This project has been written for mutual radio occultation of the Martian atmosphere. Where  radio signals are  sent from MarsExpress (MEx) to ExoMars' Trace Gas Orbiter (TGO). This can be treated as an example, where the satelites and even the planet are up to the end user's discretion. Fortunately, with the SPICE toolset, everything in the analysis can be changed with simple alterations of the starting variables.*


## Prerequisites
+ Have installed [Cosmographia 4.0](https://naif.jpl.nasa.gov/naif/cosmographia.html)
+ Some knowledge of [SPICE](https://naif.jpl.nasa.gov/naif/tutorials.html). This may take some time to understand but is very much worth it. SPICE is a remarkably powerful free-to-use tool.
+ Install the correct SPICE kernels, these are the datasets that the SPICE application reads. They include infomation on spacecraft, groundstation and planetary ephemarides, planetary shapes, onboard instrument details, leapsecond tables, reference frame conversions and *much more*. Depending on the mission the user is interested in, the download location will vary. Space agencies will manage the SPICE kernels for thier own mission. [ESA](https://www.cosmos.esa.int/web/spice) manage the kernels for [MEx](https://www.cosmos.esa.int/web/spice/spice-for-mex) and [TGO](ftp://spiftp.esac.esa.int/data/SPICE/ExoMars2016/).


## What is Radio Occultation?
In perticular this project will focus on mutual occultation, though the code could be easily edited to include conventional spacecraft->planatary atmosphere->earth occultations also. Mutual occultation (sometimes refered to as Cross-Link Occultatation) is a method of passing a  radiowave through an atmosphere between two spacecraft, to garner atmospheric parameters from the doppler shift exibted on the wave. This way the SNR is far better because of the proximity of the two spacecraft and due to there being no dispersive space and earth atmosphere in the radiowave's path. 

## Producing a ```catalog file``` for Cosmographia
To produce and image similar to the one shown at the top of the page:
1. Achieve the prerequisites
2. Change the variables to the spacecraft, planet and date range of choice:
```
#-----------------------------------------------------<VALUES TO EDIT REGULARLY>----------------------------------------
# If you only wish to analysis mutual [cross-link] occultation between MEX and TGO, then this is the only section that
# needs to be edited
start  = '2020 JAN 1'
stop   = '2020 JAN 6'
OCCSELECTION = 2 # Which occultation do you wish to see?
PathtoMetaKernel1 = '_____'
PathtoMetaKernel2 = '_____' #[ ommit if only one mk]
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
```
3. This will create a Cosmographia readable .json file, now go and ``` Open Catalog... ``` and select this new Profile.json


## Feedback
Feel free file an issue. Contributors and feature requests are always welcome. 
