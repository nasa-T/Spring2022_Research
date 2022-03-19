import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.constants import R_sun, b_wien
import astropy.units as u
import matplotlib as mpl

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

# Reading in non-tight-binary WISE data
from astropy.io import ascii
FOLDER = '/Users/Celloman/Research/Spring2022_Research/infrareddata/'
NTFILE1 = FOLDER+'DaRio16_wise.txt'
NTFILE2 = FOLDER+'Kounkel16_wise.txt'

DaRwise = ascii.read(NTFILE1)
Kwise = ascii.read(NTFILE2)
# Data with high quality:
DaRwHQ = DaRwise[DaRwise['ph_qual']=='AAAA']
KwHQ = Kwise[Kwise['ph_qual']=='AAAA']

W1DaR = DaRwHQ['w1mpro']
W2DaR = DaRwHQ['w2mpro']
W3DaR = DaRwHQ['w3mpro']
W4DaR = DaRwHQ['w4mpro']
W1K = KwHQ['w1mpro']
W2K = KwHQ['w2mpro']
W3K = KwHQ['w3mpro']
W4K = KwHQ['w4mpro']

# Finding Average Teff of Da Rio 2016
DaRio = pd.read_csv('~/Research/Spring2022_Research/DaRio2016_Teff.csv')
DaRioTemp = np.mean(DaRio.loc[:,'Teff']) * u.K
Kounkel = pd.read_csv('~/Research/Spring2022_Research/Kounkel2016_Teff.csv')
KounkelTemp = np.mean(Kounkel.loc[:,'Teff']) * u.K
dDaRW1 = d(DaRioTemp,Tdust(W1_wl))
dDaRW2 = d(DaRioTemp,Tdust(W2_wl))
dDaRW3 = d(DaRioTemp,Tdust(W3_wl))
dDaRW4 = d(DaRioTemp,Tdust(W4_wl))
dKW1 = d(DaRioTemp,Tdust(W1_wl))
dKW2 = d(DaRioTemp,Tdust(W2_wl))
dKW3 = d(DaRioTemp,Tdust(W3_wl))
dKW4 = d(DaRioTemp,Tdust(W4_wl))

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

norm = mpl.colors.Normalize(vmin=np.min(Tstar2), vmax=np.max(Tstar2))
cmap = mpl.cm.autumn
sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
for i in range(len(df2)):
    distances = np.array([d2W1.loc[i], d2W2.loc[i], d2W3.loc[i], d2W4.loc[i]])
    mags = np.array([W1_2.loc[i], W2_2.loc[i], W3_2.loc[i], W4_2.loc[i]])
    plt.plot(distances, mags, c=sm.to_rgba(Tstar2.loc[i]))
for i in range(len(DaRwHQ)):
    distances = np.array([dDaRW1.value, dDaRW2.value, dDaRW3.value, dDaRW4.value])
    mags = np.array([W1DaR[i], W2DaR[i], W3DaR[i], W4DaR[i]])
    plt.plot(distances, mags, 'g--', alpha=0.03)
for i in range(len(KwHQ)):
    distances = np.array([dKW1.value, dKW2.value, dKW3.value, dKW4.value])
    mags = np.array([W1K[i], W2K[i], W3K[i], W4K[i]])
    plt.plot(distances, mags, 'g--', alpha=0.03)
plt.colorbar(sm, label='Teff', orientation='vertical')
plt.gca().invert_yaxis()
plt.xlabel('Distance of Emission (au)')
plt.ylabel('Magnitude of Emission')
plt.show()

