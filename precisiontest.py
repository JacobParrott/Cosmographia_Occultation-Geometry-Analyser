from mpmath import *
from decimal import *

division = 72/7
print(division)

getcontext().prec = 90 ; mp.dps = 90
divisiondecimal = Decimal(72) / Decimal(7) #cause you have imported everything you dont need the the dot notation
print(divisiondecimal)

#practice with high precsion trig

mp.pretty = True
a = mp.sin(divisiondecimal)
print(a)

print('Test complete')