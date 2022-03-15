import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.constants import R_sun, b_wien
import astropy.units as u

SOURCE = '~/Research/Spring2022_Research/ALMA_tight_binaries_Data_keepers2.csv'

df = pd.read_csv(SOURCE)
df = df.loc[:108] # non-data rows removed
W1_wl = 3.4 * u.micron
W2_wl = 4.6 * u.micron
W3_wl = 12 * u.micron
W4_wl = 22 * u.micron
Tdust = lambda wl: b_wien.cgs / wl.cgs
Tstar = df.loc[:,'Teff'] * u.K
# stars assumed to have radius ~ 2 * Rsol:
d = lambda Ts, Td: (2 * R_sun.to(u.au)) / 2 * ((Ts/Td)**2)

dW1 = d(Tstar,Tdust(W1_wl))
dW2 = d(Tstar,Tdust(W2_wl))
dW3 = d(Tstar,Tdust(W3_wl))
dW4 = d(Tstar,Tdust(W4_wl))

plt.plot(Tstar, dW1, 'o')
plt.plot(Tstar, dW2,'o')
plt.plot(Tstar, dW3,'o')
plt.plot(Tstar, dW4,'o')
plt.legend(['W1','W2','W3','W4'])
plt.xlabel('Stellar Temperature')
plt.ylabel('Distance of Emission')
plt.show()
