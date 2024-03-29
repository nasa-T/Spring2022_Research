# flux to mass: M = (F * dist**2)/ (2.3 * 0.01 * B(20))
# dust opacity = 2.3; dust-to-gas ratio = 0.01; 
# B(20) = Planck blackbody function at avg-temp = 20 K;
# wavelength = 1.3mm

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.constants import c, h, k_B, M_jup, pc

SOURCE = '~/Research/Spring2022_Research/ALMA_tight_binaries_Data_keepers2.csv'
df = pd.read_csv(SOURCE)
df = df.loc[:108] # non-data rows removed
detected = df.loc[:,'Disk Mass err'] != -10 # boolean mask for detections
# change 'null' error values to half of parallax range of Orion (2 < p < 3.5):
df['p_error'] = df['p_error'].mask(df.loc[:,'p_error'] == -10, 0.75) 
flux = df.loc[:,'Flux (Jy)'] * 1E-23
flux_err = df.loc[:,'Flux_err'] * 1E-23 # flux error
dist = 1000/df.loc[:,'p'] * pc.cgs.value # distance
dist_err = df.loc[:,'p_error'] / df.loc[:,'p'] * dist # distance error


wl = .13 # wavelength of interest in centimeters (ALMA Band 6)
freq = c.cgs.value / wl # frequency of interest
freq = (249.5 + 252 + 265 + 267.56) / 4 * 10**9
# freq = 255.5E9
z = 0.01 # dust-to-gas ratio
# opacity = 2.3 # dust opacity in cm^2/g
kappa = lambda frequency: 10 * (frequency/1E12)
opacity = kappa(freq)
avgTemp = 20 # average disk temperature in Kelvin (Andrews et al., 2013)
Mjup = M_jup.cgs.value # mass of jupiter in g

# Planck blackbody function:
B = lambda frequency, temp: (2 * h.cgs.value * (frequency ** 3)) / (c.cgs.value**2 * (np.e**(h.cgs.value*frequency/(k_B.cgs.value*temp)) - 1))
    # (2 * h.cgs.value * (c.cgs.value**2) / (wavelength**5)) / (np.e**(h.cgs.value*c.cgs.value/(wavelength*k_B.cgs.value*temp))) 
# flux to mass (in Mjup) conversion:
mass = lambda F, distance: ((F * distance**2) / (opacity * z * B(freq,20))) / Mjup
print(opacity, B(freq, 20))
# print(mass(flux, dist))
from lifelines import KaplanMeierFitter
n = 10 # how many different plots
ax = plt.subplot(1,1,1)
for i in range(n):
    # fluxes and distances randomly picked in gaussian distribution:
    rFlux = flux + flux_err * np.random.randn() 
    rDist = dist + dist_err * np.random.randn()
    masses = mass(rFlux, rDist) # masses calculated from random fluxes and distances
    # fit lifelines
    kmf = KaplanMeierFitter()
    kmf.fit(masses, detected)
    kmf.plot_cumulative_density(legend=False, ax=ax)
    ax.set_xlabel('Mass (Mjup)')
# plt.savefig('distribution_of_lifelines.png')
plt.show()
