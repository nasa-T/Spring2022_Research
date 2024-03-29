import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.constants import R_sun, b_wien, c
import astropy.units as u
import matplotlib as mpl

SOURCE = '~/Research/Spring2022_Research/ALMA_tight_binaries_Data_keepers2.csv'

df = pd.read_csv(SOURCE)
df = df.loc[:108] # non-data rows removed
J_wl = 1.235 * u.micron
H_wl = 1.662 * u.micron
K_wl = 2.159 * u.micron
W1_wl = 3.3526 * u.micron
W2_wl = 4.6028 * u.micron
W3_wl = 11.5608 * u.micron
W4_wl = 22.0883 * u.micron
# wavelength of each band
wls = np.array([ 1.235, 1.662, 2.159, 3.3526, 4.6028, 11.5608, 22.0883 ]) * u.micron
J_freq = c.cgs.to(u.micron/u.s) / J_wl
H_freq = c.cgs.to(u.micron/u.s) / H_wl
K_freq = c.cgs.to(u.micron/u.s) / K_wl
W1_freq = c.cgs.to(u.micron/u.s) / W1_wl
W2_freq = c.cgs.to(u.micron/u.s) / W2_wl
W3_freq = c.cgs.to(u.micron/u.s) / W3_wl
W4_freq = c.cgs.to(u.micron/u.s) / W4_wl
# frequency of each band
freqs = c.cgs.to(u.micron/u.s) / wls
# dust temperature using Wien's law
Tdust = lambda wl: b_wien.cgs / wl.cgs
Tstar = df.loc[:,'Teff'] * u.K
W1 = df.loc[:,'W1']
W2 = df.loc[:,'W2']
W3 = df.loc[:,'W3']
W4 = df.loc[:,'W4']
# magnitudes of bands in tight binary data
tMags = np.array([ df.loc[:,'J'], df.loc[:,'H'], df.loc[:,'K'], df.loc[:,'W1'], 
                df.loc[:,'W2'], df.loc[:,'W3'], df.loc[:,'W4'] ])
# stars assumed to have radius ~ 2 * Rsol:
d = lambda Ts, Td: (2 * R_sun.to(u.au)) / 2 * ((Ts/Td)**2)

dJ = d(Tstar,Tdust(J_wl))
dH = d(Tstar,Tdust(H_wl))
dK = d(Tstar,Tdust(K_wl))
dW1 = d(Tstar,Tdust(W1_wl))
dW2 = d(Tstar,Tdust(W2_wl))
dW3 = d(Tstar,Tdust(W3_wl))
dW4 = d(Tstar,Tdust(W4_wl))

dists = np.empty(7) * u.AU # initialize 2D array for distances
for i in range(len(df)):
    # add rows to array: each row is an object with distances for each wl
    dists = np.vstack( [dists, d(Tstar[i], Tdust(wls))] ) 
dists = dists[1:] # delete initializer row

for i in range(7):
    plt.plot(Tstar, dists[:,i],'o')
# # plt.plot(Tstar, dJ, 'o')
# # plt.plot(Tstar, dH, 'o')
# # plt.plot(Tstar, dK, 'o')
# # plt.plot(Tstar, dW1, 'o')
# # plt.plot(Tstar, dW2,'o')
# # plt.plot(Tstar, dW3,'o')
# # plt.plot(Tstar, dW4,'o')
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
DaR = join(DaR2mHQ, DaRwHQ, keys=['ra_01','dec_01'])
K = join(K2mHQ, KwHQ, keys=['ra_01','dec_01'])

# Magnitudes of each band and conversion to fluxes for DaRio and Kounkel data
JDaR = DaR['j_m']
HDaR = DaR['h_m']
KDaR = DaR['k_m']
W1DaR = DaR['w1mpro']
W2DaR = DaR['w2mpro']
W3DaR = DaR['w3mpro']
W4DaR = DaR['w4mpro']
DaRMags = np.array([ np.median(DaR['j_m']),np.median(DaR['h_m']),np.median(DaR['k_m']),np.median(DaR['w1mpro']),np.median(DaR['w2mpro']),np.median(DaR['w3mpro']),np.median(DaR['w4mpro']) ])
DaRFluxes = Fν(Fν0s, DaRMags)
JK = K['j_m']
HK = K['h_m']
KK = K['k_m']
W1K = K['w1mpro']
W2K = K['w2mpro']
W3K = K['w3mpro']
W4K = K['w4mpro']
KMags = np.array([ np.median(K['j_m']),np.median(K['h_m']),np.median(K['k_m']),np.median(K['w1mpro']),np.median(K['w2mpro']),np.median(K['w3mpro']),np.median(K['w4mpro']) ])
KFluxes = Fν(Fν0s, KMags)

# Finding Average Teff of Da Rio 2016
DaRio = pd.read_csv('~/Research/Spring2022_Research/DaRio2016_Teff.csv')
DaRioTemp = np.mean(DaRio.loc[:,'Teff']) * u.K
Kounkel = pd.read_csv('~/Research/Spring2022_Research/Kounkel2016_Teff.csv')
KounkelTemp = np.mean(Kounkel.loc[:,'Teff']) * u.K
# distances for corresponding wavelengths in DaRio and Kounkel data
dDaRJ = d(DaRioTemp,Tdust(J_wl))
dDaRH = d(DaRioTemp,Tdust(H_wl))
dDaRK = d(DaRioTemp,Tdust(K_wl))
dDaRW1 = d(DaRioTemp,Tdust(W1_wl))
dDaRW2 = d(DaRioTemp,Tdust(W2_wl))
dDaRW3 = d(DaRioTemp,Tdust(W3_wl))
dDaRW4 = d(DaRioTemp,Tdust(W4_wl))
DaRDists = d(DaRioTemp, Tdust(wls))
dKJ = d(KounkelTemp,Tdust(J_wl))
dKH = d(KounkelTemp,Tdust(H_wl))
dKK = d(KounkelTemp,Tdust(K_wl))
dKW1 = d(KounkelTemp,Tdust(W1_wl))
dKW2 = d(KounkelTemp,Tdust(W2_wl))
dKW3 = d(KounkelTemp,Tdust(W3_wl))
dKW4 = d(KounkelTemp,Tdust(W4_wl))
KDists = d(KounkelTemp, Tdust(wls))

# Making "profiles" of disks with all WISE bands
# gather only objects with data from all WISE bands:
# zero-magnitude fluxes for all bands
FνJ = 1594
FνH = 1024
FνK = 666.7
FνW1 = 306.682
FνW2 = 170.663
FνW3 = 29.045
FνW4 = 8.284
FλW1 = 8.1787E-15
FλW2 = 2.4150E-15
FλW3 = 6.5151E-17
FλW4 = 5.0901E-18
df2 = df[df.loc[:,'W1'].notnull() & df.loc[:,'W2'].notnull() & df.loc[:,'W3'].notnull() & df.loc[:,'W4'].notnull()].reset_index()
# magnitudes of bands in tight binary data
tMags2 = np.array([ df2.loc[:,'J'], df2.loc[:,'H'], df2.loc[:,'K'], 
                    df2.loc[:,'W1'], df2.loc[:,'W2'], df2.loc[:,'W3'], 
                    df2.loc[:,'W4'] ])
# tMags2 = np.empty(7) * u.AU # initialize 2D array for distances
tMags2 = np.swapaxes(tMags2,0,1)
# print(tMags2[:][:2])
tFluxes = np.empty(7)
for i in range(len(df2)):
    # add rows to array: each row is an object with fluxes for each wl
    tFluxes = np.vstack( [tFluxes, Fν(Fν0s, tMags2[:][i])] ) 
tFluxes = tFluxes[1:] # delete initializer row
Tstar2 = df2.loc[:,'Teff'] * u.K
J_2Fν = Fν(FνJ,df2.loc[:,'J'])
H_2Fν = Fν(FνH,df2.loc[:,'H'])
K_2Fν = Fν(FνK,df2.loc[:,'K'])
W1_2Fν = Fν(FνW1,df2.loc[:,'W1'])
W2_2Fν = Fν(FνW2,df2.loc[:,'W2'])
W3_2Fν = Fν(FνW3,df2.loc[:,'W3'])
W4_2Fν = Fν(FνW4,df2.loc[:,'W4'])
W1_2Fλ = Fλ(FλW1,df2.loc[:,'W1'])
W2_2Fλ = Fλ(FλW2,df2.loc[:,'W2'])
W3_2Fλ = Fλ(FλW3,df2.loc[:,'W3'])
W4_2Fλ = Fλ(FλW4,df2.loc[:,'W4'])
dJ = d(Tstar2,Tdust(J_wl))
dH = d(Tstar2,Tdust(H_wl))
dK = d(Tstar2,Tdust(K_wl))
d2W1 = d(Tstar2,Tdust(W1_wl))
d2W2 = d(Tstar2,Tdust(W2_wl))
d2W3 = d(Tstar2,Tdust(W3_wl))
d2W4 = d(Tstar2,Tdust(W4_wl))
dists2 = np.empty(7) * u.AU # initialize 2D array for distances
for i in range(len(df2)):
    # add rows to array: each row is an object with with distances for each wl
    dists2 = np.vstack( [dists2, d(Tstar2[i], Tdust(wls))] ) 
dists2 = dists2[1:] # delete initializer row

norm = mpl.colors.Normalize(vmin=np.min(Tstar2), vmax=np.max(Tstar2))
cmap = mpl.cm.autumn
sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
plt.xscale('log')
# plt.yscale('log')

# for i in range(len(df2)):
#     distances = np.array([dJ.loc[i], dH.loc[i], dK.loc[i], d2W1.loc[i], d2W2.loc[i], d2W3.loc[i], d2W4.loc[i]])
#     fluxes = np.log10(np.array([J_freq.value*J_2Fν.loc[i], H_freq.value*H_2Fν.loc[i], K_freq.value*K_2Fν.loc[i], W1_freq.value*W1_2Fν.loc[i], W2_freq.value*W2_2Fν.loc[i],W3_freq.value*W3_2Fν.loc[i], W4_freq.value*W4_2Fν.loc[i]])/(J_freq.value*J_2Fν.loc[i]))
#     plt.text(distances[-1],fluxes[-1],df2.loc[i,'Disk Mass (Mjup)'])
#     if df2.loc[i,'Disk Mass err'] != -10:
#         plt.plot(distances, fluxes, c=sm.to_rgba(Tstar2.loc[i]), label='Detection')
#     else:
#         plt.plot(distances, fluxes, '--', c=sm.to_rgba(Tstar2.loc[i]), label='Non-Detection')
for i in range(len(df2)):
    fluxes = np.log10(freqs.value*tFluxes[i]/(freqs[0].value*tFluxes[i][0]))
    # print(dists[i][-1], fluxes[-1])
    plt.text(dists2[i][-1].value,fluxes[-1],df2.loc[i,'Disk Mass (Mjup)'])
    if df2.loc[i,'Disk Mass err'] != -10:
        plt.plot(dists2[i], fluxes, c=sm.to_rgba(Tstar2.loc[i]), label='Detection')
    else:
        plt.plot(dists2[i], fluxes, '--', c=sm.to_rgba(Tstar2.loc[i]), label='Non-Detection')

distancesDaR = DaRDists.value
# np.array([dDaRJ.value, dDaRH.value, dDaRK.value, dDaRW1.value, dDaRW2.value, dDaRW3.value, dDaRW4.value])
fluxesDaR = np.log10(freqs*DaRFluxes/(freqs[0]*DaRFluxes[0]))
# np.log10(np.array([J_freq.value*Fν(FνJ,np.median(JDaR)), H_freq.value*Fν(FνH,np.median(HDaR)), K_freq.value*Fν(FνK,np.median(KDaR)), W1_freq.value*Fν(FνW1,np.median(W1DaR)), W2_freq.value*Fν(FνW2,np.median(W2DaR)), W3_freq.value*Fν(FνW3,np.median(W3DaR)), W4_freq.value*Fν(FνW4,np.median(W4DaR))])/(J_freq.value*Fν(FνJ,np.median(JDaR))))
fDaRstDev = np.sqrt(np.abs(fluxesDaR))
plt.plot(distancesDaR, fluxesDaR, 'g--', alpha=0.4, linewidth=5)
plt.fill_between(distancesDaR,fluxesDaR-fDaRstDev,fluxesDaR+fDaRstDev,alpha=0.2)

distancesK = KDists.value
# np.array([dKJ.value, dKH.value, dKK.value, dKW1.value, dKW2.value, dKW3.value, dKW4.value])
fluxesK = np.log10(freqs*KFluxes/(freqs[0]*KFluxes[0]))
# np.log10(np.array([J_freq.value*Fν(FνJ,np.median(JK)), H_freq.value*Fν(FνH,np.median(HK)), K_freq.value*Fν(FνK,np.median(KK)), W1_freq.value*Fν(FνW1,np.median(W1K)), W2_freq.value*Fν(FνW2,np.median(W2K)), W3_freq.value*Fν(FνW3,np.median(W3K)), W4_freq.value*Fν(FνW4,np.median(W4K))])/(J_freq.value*Fν(FνJ,np.median(JK))))
fKstDev = np.sqrt(np.abs(fluxesK))
plt.plot(distancesK, fluxesK, 'b--', alpha=0.4, linewidth=5)
plt.fill_between(distancesK,fluxesK-fKstDev,fluxesK+fKstDev,alpha=0.2)

# for i in range(len(DaR)):
#     distances = np.array([dDaRJ.value, dDaRH.value, dDaRK.value, dDaRW1.value, dDaRW2.value, dDaRW3.value, dDaRW4.value])
#     mags = np.log10(np.array([J_freq.value*Fν(FνJ,JDaR[i]), H_freq.value*Fν(FνH,HDaR[i]), K_freq.value*Fν(FνK,KDaR[i]), W1_freq.value*Fν(FνW1,W1DaR[i]), W2_freq.value*Fν(FνW2,W2DaR[i]), W3_freq.value*Fν(FνW3,W3DaR[i]), W4_freq.value*Fν(FνW4,W4DaR[i])])/(J_freq.value*Fν(FνJ,JDaR[i])))
#     plt.plot(distances, mags, 'g--', alpha=0.03)

# for i in range(len(K)):
#     distances = np.array([dKJ.value, dKH.value, dKK.value, dKW1.value, dKW2.value, dKW3.value, dKW4.value])
#     mags = np.log10(np.array([J_freq.value*Fν(FνJ,JK[i]), H_freq.value*Fν(FνH,HK[i]), K_freq.value*Fν(FνK,KK[i]), W1_freq.value*Fν(FνW1,W1K[i]), W2_freq.value*Fν(FνW2,W2K[i]), W3_freq.value*Fν(FνW3,W3K[i]), W4_freq.value*Fν(FνW4,W4K[i])])/(J_freq.value*Fν(FνJ,JK[i])))
#     plt.plot(distances, mags, 'g--', alpha=0.03)

plt.colorbar(sm, label='Teff', orientation='vertical')
# plt.gca().invert_yaxis()
plt.title('Light Curve of Disks')
plt.xlabel('Distance of Emission (au)')
plt.ylabel('Flux of Emission log( (νFν) / (νFν)J )')

import matplotlib.lines as mlines
legendLine1 = mlines.Line2D([], [], color='black', linestyle='-', label='Detections')
legendLine2 = mlines.Line2D([], [], color='black', linestyle='--', label='Non-Detections')
# print(legendLine1)
# plt.legend(handles=legendLine1)
plt.legend([legendLine1,legendLine2],['Detections','Non-Detections'])
plt.show()


