import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SOURCE = '~/Research/Spring2022_Research/ALMA_tight_binaries_Data_keepers2.csv'

df = pd.read_csv(SOURCE)
df = df.loc[:108] # non-data rows removed
# JHK data for tight-binaries
Jt = df.loc[:,'J']
Ht = df.loc[:,'H']
Kt = df.loc[:,'K']
W1t = df.loc[:,'W1']
W2t = df.loc[:,'W2']
W3t = df.loc[:,'W3']
W4t = df.loc[:,'W4']

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
KwHQ = Kwise[Kwise['ph_qual']=='AAAA']
# Plot non-tight-binary data:
plt.plot(DaR2mHQ['h_k'], DaR2mHQ['j_h'], 'ro', markersize=0.5)
plt.plot(K2mHQ['h_k'], K2mHQ['j_h'], 'ko', markersize=0.5)
# CTTS locus equation:
j_h = lambda h, k: 0.58 * (h-k) + 0.52

# Plot tight-binary data:
plt.plot(Ht-Kt, Jt-Ht, 'o', markersize=4)
# Plot CTTS locus:
plt.plot(Ht-Kt, j_h(Ht,Kt), 'k-')
# Interstellar reddening vectors:
ISR_vector = lambda x, c: 1.83 * x + c
x1 = np.linspace(1.02/1.25,2.25,2)
x2 = np.linspace(0.618/1.25,2.25,2)
# Make dataframe for objects above ISR vector; Class III YSOs:
JHKclassIIIDaR = DaR2mHQ[DaR2mHQ['j_h'] > ISR_vector(DaR2mHQ['h_k'],-0.098)]
JHKclassIIIK = K2mHQ[K2mHQ['j_h'] > ISR_vector(K2mHQ['h_k'],-0.098)]
JHKclassIIIt = df[(Jt-Ht) > ISR_vector(Ht-Kt,-0.098)]
# Make dataframe for objects below ISR vector; Class II YSOs:
JHKclassIIDaR = DaR2mHQ[DaR2mHQ['j_h'] < ISR_vector(DaR2mHQ['h_k'],-0.098)]
JHKclassIIK = K2mHQ[K2mHQ['j_h'] < ISR_vector(K2mHQ['h_k'],-0.098)]
JHKclassIIt = df[(Jt-Ht) < ISR_vector(Ht-Kt,-0.098)]
print("Class II Percentage for tight binary sample (JHK):",len(JHKclassIIt)/len(JHKclassIIIt) * 100, "%")
print("For DaRio 2016 sample:", len(JHKclassIIDaR)/len(JHKclassIIIDaR) * 100, "%")
print("For Kounkel 2016 sample:", len(JHKclassIIK)/len(JHKclassIIIK) * 100, "%")
# c values found from Robberto et al. (2010):
# plt.plot(x1, ISR_vector(x1, -0.5),'k-') 
plt.plot(x2, ISR_vector(x2, -0.098),'k-')
plt.legend(['DaRio2016', 'Kounkel2016', 'Tight-Binaries'])
plt.show()
W23DaR = DaRwHQ['w2mpro']-DaRwHQ['w3mpro']
W12DaR = DaRwHQ['w1mpro']-DaRwHQ['w2mpro']
W23K = KwHQ['w2mpro']-KwHQ['w3mpro']
W12K = KwHQ['w1mpro']-KwHQ['w2mpro']
W23t = W2t-W3t
W12t = W1t-W2t
plt.plot(W23DaR, W12DaR,'o', markersize=0.5)
plt.plot(W23K, W12K,'o', markersize=0.5)
plt.plot(W23t, W12t, 'o', markersize=4)
plt.legend(['DaRio2016', 'Kounkel2016', 'Tight-Binaries'])

# Boundaries for YSOs
y3 = lambda Δw: 0.46 * (Δw) - 0.9 # also eqn 8
y4 = lambda Δw: -0.42 * (Δw) + 2.2 # also eqn 9
y6 = lambda Δw: 0.71 * (Δw) - 0.07
y7 = lambda Δw: -1.5 * (Δw) + 2.1

plt.plot(np.ones(10)*2.0, np.linspace(1.3,3.5,10), 'k-')
plt.plot(np.ones(10)*4.5, np.linspace(1.15,3.5,10), 'k-')
plt.plot(np.linspace(2.5,4.5,10), y3(np.linspace(2.5,4.5,10)), 'k-')
plt.plot(np.linspace(2,3.5,10), y4(np.linspace(2,3.5,10)), 'k-')

plt.plot(np.linspace(1.25,2.5,10), np.ones(10)*0.25, 'k-')
plt.plot(np.linspace(0.99,2,10),y6(np.linspace(0.99,2,10)), 'k-')
plt.plot(np.linspace(0.99,1.25,10), y7(np.linspace(0.99,1.25,10)),'k-')

WclassIIt = df[ (W12t > 0.25) & (W12t < y6(W23t)) & (W12t > y7(W23t)) & (W12t > y3(W23t)) & (W12t < y4(W23t)) ]
WclassIIDaR = DaRwHQ[ (W12DaR > 0.25) & (W12DaR < y6(W23DaR)) & (W12DaR > y7(W23DaR)) & (W12DaR > y3(W23DaR)) & (W12DaR < y4(W23DaR)) ]
WclassIIK = KwHQ[ (W12K > 0.25) & (W12K < y6(W23K)) & (W12K > y7(W23K)) & (W12K > y3(W23K)) & (W12K < y4(W23K)) ]

print("Class II percentage of sample for tight binary sample (WISE 1,2,3):",len(WclassIIt)/len(df[W1t.notnull() & W2t.notnull() & W3t.notnull()]) * 100, "%")
print("For DaRio 2016 sample:", len(WclassIIDaR)/len(DaRwHQ) * 100, "%")
print("For Kounkel 2016 sample:", len(WclassIIK)/len(KwHQ) * 100, "%")
plt.show()
