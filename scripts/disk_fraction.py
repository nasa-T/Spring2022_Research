import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SOURCE = '~/Research/Spring2022_Research/ALMA_tight_binaries_Data_keepers2.csv'

df = pd.read_csv(SOURCE)
df = df.loc[:108] # non-data rows removed
df2 = df[df['Disk Mass err']!=-10]
df3 = df[df['Disk Mass err']==-10]
# JHK+WISE data for tight-binaries
Jt = df.loc[:,'J']
Ht = df.loc[:,'H']
Kt = df.loc[:,'K']
W1t = df.loc[:,'W1']
W2t = df.loc[:,'W2']
W3t = df.loc[:,'W3']
W4t = df.loc[:,'W4']
Jt2 = df2.loc[:,'J']
Ht2 = df2.loc[:,'H']
Kt2 = df2.loc[:,'K']
W1t2 = df2.loc[:,'W1']
W2t2 = df2.loc[:,'W2']
W3t2 = df2.loc[:,'W3']
W4t2 = df2.loc[:,'W4']
Jt3 = df3.loc[:,'J']
Ht3 = df3.loc[:,'H']
Kt3 = df3.loc[:,'K']
W1t3 = df3.loc[:,'W1']
W2t3 = df3.loc[:,'W2']
W3t3 = df3.loc[:,'W3']
W4t3 = df3.loc[:,'W4']
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

# W3-W4, W2-W3, W1-W2
W34DaR = DaRwHQ['w3mpro']-DaRwHQ['w4mpro']
W23DaR = DaRwHQ['w2mpro']-DaRwHQ['w3mpro']
W12DaR = DaRwHQ['w1mpro']-DaRwHQ['w2mpro']
W34K = KwHQ['w3mpro']-KwHQ['w4mpro']
W23K = KwHQ['w2mpro']-KwHQ['w3mpro']
W12K = KwHQ['w1mpro']-KwHQ['w2mpro']
W34t = W3t-W4t
W23t = W2t-W3t
W12t = W1t-W2t
W34t2 = W3t2-W4t2
W23t2 = W2t2-W3t2
W12t2 = W1t2-W2t2
W34t3 = W3t3-W4t3
W23t3 = W2t3-W3t3
W12t3 = W1t3-W2t3

# CTTS locus equation:
j_h = lambda h, k: 0.58 * (h-k) + 0.52
# MS Curve from Pecaut & Mamajek 2013:
MSjh = np.array([0.11,0.14,0.15,0.16,0.17,0.19,0.21,0.22,0.23,0.23,0.24,0.25,0.27,0.29,0.31,0.33,0.35,0.37,0.39,0.41,0.43,0.46,0.49,0.55,0.60,0.64,0.66,0.66,0.67,0.67,0.68,0.68,0.67,0.66,0.62,0.55,0.54,0.58,0.65,0.70])
MShk = np.array([0.06,0.07,0.08,0.08,0.08,0.08,0.09,0.09,0.09,0.09,0.09,0.10,0.10,0.10,0.11,0.11,0.12,0.12,0.12,0.13,0.13,0.14,0.14,0.16,0.17,0.18,0.19,0.19,0.20,0.20,0.21,0.22,0.23,0.25,0.27,0.31,0.36,0.41,0.45,0.47])

# Interstellar reddening vectors:
ISR_vector = lambda x, c: 1.83 * x + c
# Make dataframe for objects above ISR vector; Class III YSOs:
JHKclassIIIDaR = DaR2mHQ[DaR2mHQ['j_h'] > ISR_vector(DaR2mHQ['h_k'],-0.098)]
JHKclassIIIK = K2mHQ[K2mHQ['j_h'] > ISR_vector(K2mHQ['h_k'],-0.098)]
JHKclassIIIt = df[(Jt-Ht) > ISR_vector(Ht-Kt,-0.098)]
# Make dataframe for objects below ISR vector; Class II YSOs:
JHKclassIIDaR = DaR2mHQ[DaR2mHQ['j_h'] < ISR_vector(DaR2mHQ['h_k'],-0.098)]
JHKclassIIK = K2mHQ[K2mHQ['j_h'] < ISR_vector(K2mHQ['h_k'],-0.098)]
JHKclassIIt = df[(Jt-Ht) < ISR_vector(Ht-Kt,-0.098)]

# Boundaries for YSOs (W1-W2 vs. W2-W3):
y3 = lambda Δw: 0.46 * (Δw) - 0.9 # also eqn 8
y4 = lambda Δw: -0.42 * (Δw) + 2.2 # also eqn 9
y6 = lambda Δw: 0.71 * (Δw) - 0.07
y7 = lambda Δw: -1.5 * (Δw) + 2.1

# Find Class II YSOs (W1-W2 vs. W2-W3):
WclassIIt = df[ (W12t > 0.25) & (W12t < y6(W23t)) & (W12t > y7(W23t)) & (W12t > y3(W23t)) & (W12t < y4(W23t)) ]
WclassIIDaR = DaRwHQ[ (W12DaR > 0.25) & (W12DaR < y6(W23DaR)) & (W12DaR > y7(W23DaR)) & (W12DaR > y3(W23DaR)) & (W12DaR < y4(W23DaR)) ]
WclassIIK = KwHQ[ (W12K > 0.25) & (W12K < y6(W23K)) & (W12K > y7(W23K)) & (W12K > y3(W23K)) & (W12K < y4(W23K)) ]

print("Class II Percentage for tight binary sample (JHK):",'%.2f'%(len(JHKclassIIt)/(len(JHKclassIIIt) + len(JHKclassIIt)) * 100), "±", 
'%.2f'%(np.sqrt(len(JHKclassIIt))/len(df[Jt.notnull() & Ht.notnull() & Kt.notnull()]) * 100), "%")
print("For DaRio 2016 sample:", '%.2f'%(len(JHKclassIIDaR)/len(DaR2mHQ) * 100), "±", '%.2f'%(np.sqrt(len(JHKclassIIDaR))/len(DaR2mHQ) * 100), "%")
print("For Kounkel 2016 sample:", '%.2f'%(len(JHKclassIIK)/len(K2mHQ) * 100), "±", '%.2f'%(np.sqrt(len(JHKclassIIK))/len(K2mHQ) * 100), "%")

print("Class II percentage of sample for tight binary sample (WISE 1,2,3):",
'%.2f'%(len(WclassIIt)/len(df[W1t.notnull() & W2t.notnull() & W3t.notnull()]) * 100), "±", 
'%.2f'%(np.sqrt(len(WclassIIt))/len(df[W1t.notnull() & W2t.notnull() & W3t.notnull()]) * 100), "%")
print("For DaRio 2016 sample:", '%.2f'%(len(WclassIIDaR)/len(DaRwHQ) * 100), 
"±", '%.2f'%(np.sqrt(len(WclassIIDaR))/len(DaRwHQ) * 100), "%")
print("For Kounkel 2016 sample:", '%.2f'%(len(WclassIIK)/len(KwHQ) * 100), "±", 
'%.2f'%(np.sqrt(len(WclassIIK))/len(KwHQ) * 100), "%")

fig, (ax1, ax2, ax3) = plt.subplots(1, 3)

# Plot non-tight-binary data:
ax1.plot(DaR2mHQ['h_k'], DaR2mHQ['j_h'], 'o', markersize=0.5)
ax1.plot(K2mHQ['h_k'], K2mHQ['j_h'], 'o', markersize=0.5)
# Plot tight-binary data:
ax1.plot(Ht2-Kt2, Jt2-Ht2, 'go', markersize=4)
ax1.plot(Ht3-Kt3, Jt3-Ht3, 'rx', markersize=4)
# Plot CTTS locus:
ax1.plot(Ht-Kt, j_h(Ht,Kt), 'k-')
# Plot MS curve:
ax1.plot(MShk, MSjh, 'k-')
# Plot reddening vectors:
x1 = np.linspace(1.02/1.25,2.25,2)
x2 = np.linspace(0.618/1.25,2.25,2)
# c values found from Robberto et al. (2010):
# plt.plot(x1, ISR_vector(x1, -0.5),'k-') 
ax1.plot(x2, ISR_vector(x2, -0.098),'k-')
ax1.legend(['DaRio2016', 'Kounkel2016', 'Detections', 'Non-Detections'])
ax1.set_xlabel('H-K')
ax1.set_ylabel('J-H')

# Plot W1-W2 vs. W2-W3:
ax2.plot(W23DaR, W12DaR,'o', markersize=0.5)
ax2.plot(W23K, W12K,'o', markersize=0.5)
ax2.plot(W23t2, W12t2, 'o', markersize=4)
ax2.plot(W23t3, W12t3, 'x', markersize=4)
ax2.legend(['DaRio2016', 'Kounkel2016', 'Detections', 'Non-Detections'])
# Plot boundaries found in Fischer et al. (2016):
ax2.plot(np.ones(10)*2.0, np.linspace(1.3,3.5,10), 'k-')
ax2.plot(np.ones(10)*4.5, np.linspace(1.15,3.5,10), 'k-')
ax2.plot(np.linspace(2.5,4.5,10), y3(np.linspace(2.5,4.5,10)), 'k-')
ax2.plot(np.linspace(2,3.5,10), y4(np.linspace(2,3.5,10)), 'k-')

ax2.plot(np.linspace(1.25,2.5,10), np.ones(10)*0.25, 'k-')
ax2.plot(np.linspace(0.99,2,10),y6(np.linspace(0.99,2,10)), 'k-')
ax2.plot(np.linspace(0.99,1.25,10), y7(np.linspace(0.99,1.25,10)),'k-')
ax2.set_xlabel('W2-W3')
ax2.set_ylabel('W1-W2')


# Plot W1-W2 vs. W3-W4:
ax3.plot(W34DaR, W12DaR,'o', markersize=0.5)
ax3.plot(W34K, W12K,'o', markersize=0.5)
ax3.plot(W34t2, W12t2, 'o', markersize=4)
ax3.plot(W34t3, W12t3, 'x', markersize=4)
ax3.legend(['DaRio2016', 'Kounkel2016', 'Detections', 'Non-Detections'])
ax3.set_xlabel('W3-W4')
ax3.set_ylabel('W1-W2')

plt.show()
