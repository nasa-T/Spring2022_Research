"""dust_distance.py...but better"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.constants import R_sun, b_wien, c
import astropy.units as u
import matplotlib as mpl

SOURCE = '~/Research/Spring2022_Research/ALMA_tight_binaries_Data_keepers2.csv'

df = pd.read_csv(SOURCE)
df = df.loc[:108] # non-data rows removed
# wavelength of each band (J,H,K,W1-4)
wls = np.array([ 1.235, 1.662, 2.159, 3.3526, 4.6028, 11.5608, 22.0883 ]) * u.micron
# frequency of each band
freqs = c.cgs.to(u.micron/u.s) / wls
# dust temperature using Wien's law
Tdust = lambda wl: b_wien.cgs / wl.cgs
Tstar = df.loc[:,'Teff'] * u.K
# magnitudes of bands in tight binary data
tMags = np.array([ df.loc[:,'J'], df.loc[:,'H'], df.loc[:,'K'], df.loc[:,'W1'], 
                df.loc[:,'W2'], df.loc[:,'W3'], df.loc[:,'W4'] ])
# stars assumed to have radius ~ 2 * Rsol:
d = lambda Ts, Td: (2 * R_sun.to(u.au)) / 2 * ((Ts/Td)**2)

dists = np.empty(7) * u.AU # initialize 2D array for distances
for i in range(len(df)):
    # add rows to array: each row is an object with distances for each wl
    dists = np.vstack( [dists, d(Tstar[i], Tdust(wls))] ) 
dists = dists[1:] # delete initializer row
# dists has dimensions Nobj x Nband (Nobj: number of objects, Nband: number of bands)

# plot distribution of distances given temps of stars
for i in range(7):
    plt.plot(Tstar, dists[:,i],'o')
plt.legend(['J','H','K','W1','W2','W3','W4'])
plt.xlabel('Stellar Temperature (K)')
plt.ylabel('Distance of Emission (au)')
plt.show()

fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2)
# plot distribution of WISE magnitudes vs. distances
ax1.scatter(dists[:,3],tMags[3])
ax2.scatter(dists[:,4],tMags[4])
ax3.scatter(dists[:,5],tMags[5])
ax4.scatter(dists[:,6],tMags[6])
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

# conversion to Fν and Fλ given zero-mag fluxes and magnitude of object in relevant band
Fν = lambda Fν0, m_vega: Fν0*10**(-m_vega/2.5)
Fλ = lambda Fλ0, m_vega: Fλ0*10**(-m_vega/2.5)
# zero-magnitude fluxes (ν) for all bands
Fν0s = np.array([1594, 1024, 666.7, 306.682, 170.663, 29.045, 8.284])

# Reading in non-tight-binary WISE data
from astropy.io import ascii
from astropy.table import join
FOLDER = '/Users/Celloman/Research/Spring2022_Research/infrareddata/'
NTFILE1 = FOLDER+'DaRio16_OrionA_2mass.txt'
NTFILE2 = FOLDER+'DaRio16_wise.txt'
NTFILE3 = FOLDER+'Kounkel16_2mass.txt'
NTFILE4 = FOLDER+'Kounkel16_wise.txt'

DaR2mass = ascii.read(NTFILE1)
DaRwise = ascii.read(NTFILE2)
K2mass = ascii.read(NTFILE3)
Kwise = ascii.read(NTFILE4)
# Data with high quality:
DaR2mHQ = DaR2mass[DaR2mass['ph_qual']=='AAA']
DaRwHQ = DaRwise[DaRwise['ph_qual']=='AAAA']
K2mHQ = K2mass[K2mass['ph_qual']=='AAA']
KwHQ = Kwise[Kwise['ph_qual']=='AAAA']
# put WISE and 2MASS data into one table for DaRio and Kounkel data
DaR = join(DaR2mHQ, DaRwHQ, keys=['ra_01','dec_01'])
K = join(K2mHQ, KwHQ, keys=['ra_01','dec_01'])

# Median magnitude of each band and conversion to fluxes for DaRio and Kounkel data
DaRMags = np.array([ np.median(DaR['j_m']),np.median(DaR['h_m']),np.median(DaR['k_m']),np.median(DaR['w1mpro']),np.median(DaR['w2mpro']),np.median(DaR['w3mpro']),np.median(DaR['w4mpro']) ])
DaRFluxes = Fν(Fν0s, DaRMags)
KMags = np.array([ np.median(K['j_m']),np.median(K['h_m']),np.median(K['k_m']),np.median(K['w1mpro']),np.median(K['w2mpro']),np.median(K['w3mpro']),np.median(K['w4mpro']) ])
KFluxes = Fν(Fν0s, KMags) # dimensions 1 x 7

# Finding Average Teff of Da Rio 2016
DaRio = pd.read_csv('~/Research/Spring2022_Research/DaRio2016_Teff.csv')
DaRioTemp = np.mean(DaRio.loc[:,'Teff']) * u.K
Kounkel = pd.read_csv('~/Research/Spring2022_Research/Kounkel2016_Teff.csv')
KounkelTemp = np.mean(Kounkel.loc[:,'Teff']) * u.K
# distances for corresponding wavelengths in DaRio and Kounkel data
# 1-D arrays
DaRDists = d(DaRioTemp, Tdust(wls))
KDists = d(KounkelTemp, Tdust(wls))

# Making "profiles" of disks with all WISE bands
# gather only objects with data from all WISE bands:
# zero-magnitude fluxes for all bands
df2 = df[df.loc[:,'W1'].notnull() & df.loc[:,'W2'].notnull() & df.loc[:,'W3'].notnull() & df.loc[:,'W4'].notnull()].reset_index()
# magnitudes of bands in tight binary data
tMags2 = np.array([ df2.loc[:,'J'], df2.loc[:,'H'], df2.loc[:,'K'], 
                    df2.loc[:,'W1'], df2.loc[:,'W2'], df2.loc[:,'W3'], 
                    df2.loc[:,'W4'] ])
tMags2 = np.swapaxes(tMags2,0,1) # make dimensions Nobj x Nband
tFluxes = np.empty(7) # initialize array for tight-binary fluxes
for i in range(len(df2)):
    # add rows to array: each row is an object with fluxes for each wl
    tFluxes = np.vstack( [tFluxes, Fν(Fν0s, tMags2[:][i])] ) 
tFluxes = tFluxes[1:] # delete initializer row
# tFluxes is a 2-D array with dimensions Nobj x Nband
Tstar2 = df2.loc[:,'Teff'] * u.K
dists2 = np.empty(7) * u.AU # initialize 2D array for distances
for i in range(len(df2)):
    # add rows to array: each row is an object with with distances for each wl
    dists2 = np.vstack( [dists2, d(Tstar2[i], Tdust(wls))] ) 
dists2 = dists2[1:] # delete initializer row

# setting up colorbar for representing temperatures of objects
norm = mpl.colors.Normalize(vmin=np.min(Tstar2), vmax=np.max(Tstar2))
cmap = mpl.cm.autumn
sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
plt.xscale('log')
# plt.yscale('log')

# loop through each object
for i in range(len(df2)):
    # make array of (normalized) fluxes for one object
    # normalization of fluxes: ν*Fν/(J-band ν*Fν)
    fluxes = np.log10(freqs.value*tFluxes[i]/(freqs[0].value*tFluxes[i][0]))
    # label object at the end of light curve with disk mass (from ALMA data)
    plt.text(dists2[i][-1].value,fluxes[-1],df2.loc[i,'Disk Mass (Mjup)'])
    # differentiate between detected and non-detected objects
    if df2.loc[i,'Disk Mass err'] != -10:
        # plot normalized fluxes vs distance for one object
        plt.plot(dists2[i], fluxes, c=sm.to_rgba(Tstar2.loc[i]), label='Detection')
    else:
        plt.plot(dists2[i], fluxes, '--', c=sm.to_rgba(Tstar2.loc[i]), label='Non-Detection')

distancesDaR = DaRDists.value # array of distances for DaR data
# normalize fluxes
fluxesDaR = np.log10(freqs*DaRFluxes/(freqs[0]*DaRFluxes[0]))
# poisson standard deviation; will change to errors provided by tables
fDaRstDev = np.sqrt(np.abs(fluxesDaR))
# plot normalized fluxes vs distance
plt.plot(distancesDaR, fluxesDaR, 'g--', alpha=0.4, linewidth=5)
# plot filled area of error for whole curve
plt.fill_between(distancesDaR,fluxesDaR-fDaRstDev,fluxesDaR+fDaRstDev,alpha=0.2)

distancesK = KDists.value
fluxesK = np.log10(freqs*KFluxes/(freqs[0]*KFluxes[0]))
fKstDev = np.sqrt(np.abs(fluxesK))
plt.plot(distancesK, fluxesK, 'b--', alpha=0.4, linewidth=5)
plt.fill_between(distancesK,fluxesK-fKstDev,fluxesK+fKstDev,alpha=0.2)

plt.colorbar(sm, label='Teff', orientation='vertical')
# plt.gca().invert_yaxis()
plt.title('Light Curve of Disks')
plt.xlabel('Distance of Emission (au)')
plt.ylabel('Flux of Emission log( (νFν) / (νFν)J )')

# create a legend for detections and non-detections using dummy lines
import matplotlib.lines as mlines
legendLine1 = mlines.Line2D([], [], color='black', linestyle='-', label='Detections')
legendLine2 = mlines.Line2D([], [], color='black', linestyle='--', label='Non-Detections')
plt.legend([legendLine1,legendLine2],['Detections','Non-Detections'])
plt.show()


