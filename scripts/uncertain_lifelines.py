# np.random.randn()
# create Lifelines based on uncertainty of flux
# convert flux to mass (Peter's Thesis)
# something trivial

# flux to mass: M = (F * dist**2)/ (2.3 * 0.01 * B(20))
# dust opacity = 2.3; dust-to-gas ratio = 0.01; 
# B(20) = Planck blackbody function at avg-temp = 20 K;
# wavelength = 1.3mm

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SOURCE = '~/Research/Spring2022_Research/ALMA_tight_binaries_Data_keepers2.csv'
df = pd.read_csv(SOURCE)
df = df.loc[:108]
flux = df.loc[:,'Flux (Jy)']
flux_err = df.loc[:,'Flux_err'] # flux error
dist = df.loc[:,'p'] # distance
dist_err = df.loc[:,'p_error'] # distance error
detected = df.loc[:,'Disk Mass err'] != -10 # boolean mask for detections
wl = 1.3E-3 # wavelength of interest (ALMA Band 6)
z = 0.01 # dust-to-gas ratio
h = 6.626E-34 # Planck's Constant
c = 2.998E8 # speed of light
kb = 1.3807E-23 # Boltzmann Constant
opacity = 2.3 # dust opacity in cm^2/g
avgTemp = 20 # average disk temperature in Kelvin (Andrews et al., 2013)
Mjup = 1.898E30 # mass of jupiter in g

# Planck blackbody function:
B = lambda wavelength, temp: (2 * h * (c**2) / (wavelength**5)) / (np.e**(h*c/(wavelength*kb*temp))) 
# flux to mass conversion:
mass = lambda F, distance: (F * distance**2) / (opacity * z * B(wl,20))

# fluxes and distances randomly picked in gaussian distribution:
rFlux = flux + flux_err * np.random.randn() 
rDist = dist + dist_err * np.random.randn() 
masses = mass(rFlux, rDist) # masses calculated from random fluxes and distances
print(masses, detected)

from lifelines import KaplanMeierFitter
kmf = KaplanMeierFitter()
kmf.fit(masses, detected)
kmf.plot_survival_function()
kmf.plot_cumulative_density()
plt.show()
