# np.random.randn()
# create Lifelines based on uncertainty of flux
# convert flux to mass (Peter's Thesis)
# something more trivial

# flux to mass: M = (F * dist**2)/ (2.3 * 0.01 * B(20))
# dust opacity = 2.3; dust-to-gas ratio = 0.01; 
# B(20) = Planck blackbody function at avg-temp = 20 K;
# wavelength = 1.3mm

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.constants import c, h, k_B, M_jup, pc
# e.g. M_jup.cgs.value
# p_errors = -10: make error range the size of cluster 

SOURCE = '~/Research/Spring2022_Research/ALMA_tight_binaries_Data_keepers2.csv'
df = pd.read_csv(SOURCE)
df = df.loc[:108] # non-data rows removed
flux = df.loc[:,'Flux (Jy)'] * 1E-23
flux_err = df.loc[:,'Flux_err'] * 1E-23 # flux error
dist = 1000/df.loc[:,'p'] * pc.cgs.value # distance
dist_err = df.loc[:,'p_error'] / df.loc[:,'p'] * dist # distance error
# print(dist_err / pc.cgs.value)

detected = df.loc[:,'Disk Mass err'] != -10 # boolean mask for detections
wl = .13 # wavelength of interest (ALMA Band 6)
freq = c.cgs.value / wl # frequency of interest
z = 0.01 # dust-to-gas ratio
opacity = 2.3 # dust opacity in cm^2/g
avgTemp = 20 # average disk temperature in Kelvin (Andrews et al., 2013)
Mjup = M_jup.cgs.value # mass of jupiter in g

# Planck blackbody function:
B = lambda frequency, temp: (2 * h.cgs.value * (frequency ** 3)) / (c.cgs.value**2 * np.e**(h.cgs.value*frequency/(k_B.cgs.value*temp)) - 1)
    # (2 * h.cgs.value * (c.cgs.value**2) / (wavelength**5)) / (np.e**(h.cgs.value*c.cgs.value/(wavelength*k_B.cgs.value*temp))) 
# flux to mass (in Mjup) conversion:
mass = lambda F, distance: ((F * distance**2) / (opacity * z * B(freq,20))) / Mjup
from lifelines import KaplanMeierFitter
n = 10 # how many different plots
for i in range(n):
    # fluxes and distances randomly picked in gaussian distribution:
    rFlux = flux + flux_err * np.random.randn() 
    rDist = dist + dist_err * np.random.randn()
    masses = mass(rFlux, rDist) # masses calculated from random fluxes and distances
    # print(masses)
    kmf = KaplanMeierFitter()
    kmf.fit(masses, detected)
    # kmf.plot_survival_function()
    kmf.plot_cumulative_density()
plt.show()
