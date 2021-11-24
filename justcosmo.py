# just run cosmo
import spiceypy as spice
from os import *
import main
A = 657605100
B = 658560507
C = 659891734
D = 660118099
et = A
here = path.abspath(path.dirname(__file__))
PathtoMetaKernel1 = here + '/TGO/mk/em16_plan.tm'
PathtoMetaKernel2 = here + '/MEX/krns/mk/MEX_OPS.tm'
spice.furnsh(PathtoMetaKernel1)
spice.furnsh(PathtoMetaKernel2)
sv = main.SpiceVariables()
main.CosmographiaCatalogFormer(et, sv)
JulianTime = spice.timout(et, sv.TFMT)
print('Gas Volumes Visualised for A ' + JulianTime)
