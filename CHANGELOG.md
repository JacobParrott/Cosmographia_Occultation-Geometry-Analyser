
# Change Log
## 26/02/20
+ Change Log added
+ ReadMe updated to reflect new capabilities
## 25/02/20

+ Hundreds of occultations can be analysed with a single run and their geometric parameters are stored into a dataframe.
+ New geometric parameters can be calculated:
	+ Solar Zenith Angle - where the sun is at the time of occultation.
	+ Grazing Angle - angle between the Martian surface's normal and the lowest -> highest point vector.
	+ Distance - Displacement between the target and the observer.
	+ Coordinates - The latitude and longitude of the lowest tangent point (at the final epoch of occultation.)
+ Multiple locations can be plotted onto a map of mars.
+ Improvements to the Cosmographia profile former. Specifically, in 20% of cases, a maximum height would never be found and the shape of the profile would become heavily distorted. Now a max height is always found.
+ Percentage of progress through analysis is printed into the users console.
+ Code is separated into a main and setup python files to conform to best coding practice.
+ Relative file paths have been added to point to the stored kernels and pictures .

## 04/02/20
+ README file added


##  29/01/20 - Initial Commit
+ Code can produce a ```catalog file``` for Cosmographia so that the profile formed from a  single mutual radio occultation can be easily visualised.
