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
W1 = df.loc[:,'W1']
W2 = df.loc[:,'W2']
W3 = df.loc[:,'W3']
W4 = df.loc[:,'W4']
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
plt.xlabel('Stellar Temperature (K)')
plt.ylabel('Distance of Emission (au)')
plt.show()

fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2)
ax1.scatter(dW1,W1)
ax2.scatter(dW2,W2)
ax3.scatter(dW3,W3)
ax4.scatter(dW4,W4)
ax1.invert_yaxis()
ax2.invert_yaxis()
ax3.invert_yaxis()
ax4.invert_yaxis()
ax1.set_title('W1')
ax2.set_title('W2')
ax3.set_title('W3')
ax4.set_title('W4')
ax1.set_xlabel('Distance of Emission (au)')
ax1.set_ylabel('Magnitude of Emission')
ax2.set_xlabel('Distance of Emission (au)')
ax2.set_ylabel('Magnitude of Emission')
ax3.set_xlabel('Distance of Emission (au)')
ax3.set_ylabel('Magnitude of Emission')
ax4.set_xlabel('Distance of Emission (au)')
ax4.set_ylabel('Magnitude of Emission')
plt.show()

# Making "profiles" of disks with all WISE bands
# gather only objects with data from all WISE bands:
df2 = df[df.loc[:,'W1'].notnull() & df.loc[:,'W2'].notnull() & df.loc[:,'W3'].notnull() & df.loc[:,'W4'].notnull()].reset_index()
Tstar2 = df2.loc[:,'Teff'] * u.K
W1_2 = df2.loc[:,'W1']
W2_2 = df2.loc[:,'W2']
W3_2 = df2.loc[:,'W3']
W4_2 = df2.loc[:,'W4']
d2W1 = d(Tstar2,Tdust(W1_wl))
d2W2 = d(Tstar2,Tdust(W2_wl))
d2W3 = d(Tstar2,Tdust(W3_wl))
d2W4 = d(Tstar2,Tdust(W4_wl))

for i in range(len(df2)):
    distances = np.array([d2W1.loc[i], d2W2.loc[i], d2W3.loc[i], d2W4.loc[i]])
    mags = np.array([W1_2.loc[i], W2_2.loc[i], W3_2.loc[i], W4_2.loc[i]])
    plt.plot(distances, mags)
plt.gca().invert_yaxis()
plt.xlabel('Distance of Emission (au)')
plt.ylabel('Magnitude of Emission')
plt.show()

