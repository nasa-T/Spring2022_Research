import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter

SOURCE = '~/Research/Spring2022_Research/ALMA_tight_binaries_Data_keepers2.csv'
df = pd.read_csv(SOURCE)
df = df.loc[:108] # non-data rows removed
detected = df.loc[:,'Disk Mass err'] != -10 # boolean mask for detections
flux = df.loc[:,'Flux (Jy)']
def rndm(a, b, g, size=1):
    """Power-law gen for pdf(x)\propto x^{g-1} for a<=x<=b
    Source: https://stackoverflow.com/questions/31114330/python-generating-random-numbers-from-a-power-law-distribution"""
    r = np.random.random(size=size)
    ag, bg = a**g, b**g
    return (ag + (bg - ag)*r)**(1./g)
# sorted minimum detection fluxes for each object (5σ)
F0 = np.sort(df.loc[:,'Flux_err']*5) 
repetitions = 9
α_range = np.linspace(-2,-0.05, repetitions) # range of α values in question
def fake_v_realSample(realFlux,realDetected,fakeFlux,fakeDetected, axis,alpha='Fake Data'):
    kmf = KaplanMeierFitter()
    kmf.fit(fakeFlux, fakeDetected, label='α = ' +str(alpha))
    kmf.plot_cumulative_density(ax=axis)
    kmf.fit(realFlux, realDetected, label='Actual Data')
    kmf.plot_cumulative_density(ax=axis)
    axis.set_xlabel('Flux (Jy)')
# ax = plt.subplot(1,1,1)
# false_sample = np.sort(rndm(min(flux), max(flux), g=-1.1+1, size=len(df)))
# false_detected = false_sample > F0
# fake_v_realSample(flux,detected,false_sample,false_detected,ax)
# plt.show()
def closestFactorPair(num):
    a = int(np.sqrt(num))
    while (num % a) != 0:
        a += 1
    return (a, num//a)
def findingAlpha(alphaRange,n=6):
    for i in range(1,n+1):
        α = round(α_range[i-1],2)
        a, b = closestFactorPair(n)
        ax = plt.subplot(a,b,i)
        false_sample = np.sort(rndm(min(flux), max(flux), g=α+1, size=len(df)))
        false_detected = false_sample > F0
        fake_v_realSample(flux,detected,false_sample,false_detected,ax, alpha=α)
    plt.show()
findingAlpha(α_range,repetitions)
# use randomized order for fake detection mask


# for i in range(1, n+1):

# kmf = KaplanMeierFitter()
# for i in range(1,n+1):
# # for α in α_range:
#     α = round(α_range[i-1],2)
#     ax = plt.subplot(3,2,i)
#     # making a plot for each α value giving a false sample (based on power-law)
#     # sort the false sample to compare to the sorted minimum detections
#     false_sample = np.sort(rndm(min(flux), max(flux), g=α+1, size=len(df)))
#     false_detected = false_sample > F0
#     kmf.fit(false_sample, false_detected, label='α = ' +str(α))
#     kmf.plot_cumulative_density(ax=ax)
#     kmf.fit(flux, detected, label='Actual Data')
#     kmf.plot_cumulative_density(ax=ax)
#     ax.set_xlabel('Flux (Jy)')
# plt.show()

