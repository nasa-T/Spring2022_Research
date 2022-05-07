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
F0 = min(df[detected].loc[:,'Flux (Jy)']) # flux of lowest detection for censoring
n = 6 # when changing this, check dimensions of subplots
α_range = np.linspace(-2,-0.05, n) # range of α values in question
kmf = KaplanMeierFitter()
for i in range(1,n+1):
# for α in α_range:
    α = round(α_range[i-1],2)
    ax = plt.subplot(3,2,i)
    # making a plot for each α value giving a false sample (based on power-law)
    false_sample = rndm(min(flux), max(flux), g=α+1, size=len(df)) 
    false_detected = false_sample > F0
    kmf.fit(false_sample, false_detected, label='α = ' +str(α))
    kmf.plot_cumulative_density(ax=ax)
    kmf.fit(flux, detected, label='Actual Data')
    kmf.plot_cumulative_density(ax=ax)
    ax.set_xlabel('Flux (Jy)')
plt.show()

