import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SOURCE = '~/Research/Spring2022_Research/ALMA_tight_binaries_Data_keepers2.csv'

df = pd.read_csv(SOURCE)
df = df.loc[:108] # non-data rows removed
J = df.loc[:,'J']
H = df.loc[:,'H']
K = df.loc[:,'K']

# Reading in non-tight-binary data
from astropy.io import ascii
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
K2wHQ = Kwise[Kwise['ph_qual']=='AAAA']
# Plot non-tight-binary data:
plt.plot(DaR2mHQ['h_k'], DaR2mHQ['j_h'], 'ro', markersize=0.5)
plt.plot(K2mHQ['h_k'], K2mHQ['j_h'], 'ko', markersize=0.5)
# CTTS locus equation:
j_h = lambda h, k: 0.58 * (h-k) + 0.52
# Make dataframe for objects above CTTS locus:
classII = (J-H) > j_h(H,K)
df2 = df[classII]
# Plot tight-binary data:
plt.plot(H-K, J-H, 'o', markersize=4)
# Plot CTTS locus:
plt.plot(H-K, j_h(H,K), 'k-')
# Interstellar reddening vectors:
ISR_vector = lambda x, c: 1.83 * x + c
x1 = np.linspace(1.02/1.25,2.25,2)
x2 = np.linspace(0.618/1.25,2.25,2)
# c values found from Robberto et al. (2010):
plt.plot(x1, ISR_vector(x1, -0.5),'k-') 
plt.plot(x2, ISR_vector(x2, -0.098),'k-')
plt.legend(['DaRio2016', 'Kounkel2016', 'Tight-Binaries'])

plt.show()

