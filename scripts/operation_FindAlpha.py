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
def fake_v_realSample(realFlux,realDetected,fakeFlux,fakeDetected, axis,alpha):
    """Use the Kaplan Meier Fitter to plot the cumulative density of two sets 
    of data on one axis; fake data is labeled with the α parameter of the power 
    law [ N(F) = A*F^-α ]"""
    kmf = KaplanMeierFitter()
    kmf.fit(fakeFlux, fakeDetected, label='α = ' +str(alpha))
    kmf.plot_cumulative_density(ax=axis)
    kmf.fit(realFlux, realDetected, label='Actual Data')
    kmf.plot_cumulative_density(ax=axis)
    axis.set_xlabel('Flux (Jy)')

def closestFactorPair(num):
    """Find the pair of factors of num that are closest together;
    implemented here for axes sizing; 
    avoid prime numbers...those are not welcome here;
    returns pair of ints which are factors of num"""
    # start with potential factor that would work if num is a perfect square
    a = int(np.sqrt(num))
    while (num % a) != 0:
        # if it is not a factor, continue to check numbers until it is
        a += 1
    # return pair of factors found
    return (a, num//a)

def findingAlpha(alphaRange):
    """Makes a bunch of plots of fake data for each alpha value and plots 
    against real data; input an array of desired alpha values to look at; 
    fruitless function that plots"""
    # n: number of alpha values, plots, etc.
    n = len(alphaRange)
    for i in range(n):
        α = alphaRange[i] # current α
        a, b = closestFactorPair(n) # get dimensions for axes
        ax = plt.subplot(a,b,i+1) # current axis
        # generate false sample based on current α
        false_sample = np.sort(rndm(min(flux), max(flux), g=α+1, size=len(df)))
        # create mask based on minimum detection limit of real sample
        false_detected = false_sample > F0
        # plot real and fake samples
        fake_v_realSample(flux,detected,false_sample,false_detected,ax, round(α,2))

findingAlpha(α_range)
# plt.savefig('fitting_alpha_to_data.png')
plt.show()
# use randomized order for fake detection mask
def randomDetectionMask(α=-1.1,n=6):
    """Similar to findingAlpha but each iteration maintains the same α but 
    randomizes the detection limit mask on each iteration"""
    # generate false sample
    false_sample = np.sort(rndm(min(flux), max(flux), g=α+1, size=len(df)))
    for i in range(1,n+1): 
        a, b = closestFactorPair(n) # get dimensions for axes
        ax = plt.subplot(a,b,i) # current axis
        if i > 1:
            # first false sample maintains order (not-randomized);
            # all others are randomized
            np.random.shuffle(false_sample)
        # create mask based on minimum detection limit of real sample
        false_detected = false_sample > F0
        # plot real and fake samples
        fake_v_realSample(flux,detected,false_sample,false_detected,ax,α)
    
randomDetectionMask(-1.2)
# plt.savefig('vary_detection_mask.png')
plt.show()
