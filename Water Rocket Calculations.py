import math
import scipy
from scipy.integrate import romberg
import numpy

P0 = 830000 #Pa, initial internal pressure
Pa = 100000 #Pa, atmospheric pressure
V0 = 0.00085 #m^3, initial volume of air
Vtot = 0.0013 #m^3, total volume
Vw0 = Vtot - V0 #m^3, initial volume of water
gamma = 1.4 #heat capacity ratio
mv = 0.280 #kg, dry mass
g = 9.81 #m s^-2, gravitational acceleration

T = 20 #°C, atmospheric temperature
hum = 0.55 #relative humidity
Pws = 2338.8 #saturation pressure

re = 0.0045 #m, exit cross-section radius
Ae = math.pi * re**2 #m^2, exit cross-section area
rho = 1300 #kg m^-3, density of fuel

rb = 0.045 #m, bottle cross-section radius
Ab = math.pi * rb**2 #m^2, bottle cross-section area
cD = 0.3 #drag coefficient
rhoExtAir = ( (28.965*Pa) / (8.314*(T+273)) ) / 1000 #kg m^2, atmospheric density
D = 0.5 * rhoExtAir * cD * Ab

vt = math.sqrt(mv * g / D) #m s^-1, terminal velocity


def f(x):
    return math.sqrt(rho/2)/Ae * (1 / math.sqrt((P0 * (V0/x)**gamma) - Pa))
tb = romberg(f, V0, Vtot) #s, burn time


VFbar = (Vtot-V0)/tb #m^3 s^-1, average volume flow
def V(t):
    return (Vtot-V0) * t / tb + V0
def m(t):
    return mv + rho * (Vtot - V(t))
def v(t):
    return 2 * Ae * (P0 * ((V0 / V(t))**gamma) - Pa) / m(t) - g
vbf = romberg(v, 0, tb) #m s^-1, velocity at tb
ybf = (vbf * tb) / 2 #m, altitude at tb

Pbf = P0 * ((V0/Vtot)**gamma) #Pa, internal pressure at tb
rhoIntAir = ( (28.965*Pbf) / (8.314*(T+273)) ) / 1000
mair = rhoIntAir * Vtot
vmax = math.sqrt( ( 2*(Pbf-Pa) ) / rhoIntAir ) * math.log( 1 + (mair/mv ) ) + vbf #m s^-1, maximum velocity

ycf = ( (vt**2) / (2*g) ) * math.log(1 + ((vmax/vt)**2) ) #m, coasting altitude

print("water burn time: " + str(tb))
print("maximum velocity: " + str(vmax))
print("highest altitude: " + str(ybf+ycf))
