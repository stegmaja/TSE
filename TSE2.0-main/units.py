'''
Unit system and constants
'''
Msol = 1.
AU = 1e-3
yr = 1.

Rsol = AU/214.95
km = 6.685e-9*AU
sec = 3.171e-8*yr
day = yr/365.
Myr = 1e6*yr
c = 63198.*AU/yr
G = 39.42*AU**3/Msol/yr**2
kms = km/sec

kL1 = .028
kL2 = .028
k1 = .1
k2 = .1
Dt = sec
vk_sigma = 265.*kms
atol, rtol = 1e-3,1e-3
step_t = 1e4*yr
tphysf = 100.*Myr