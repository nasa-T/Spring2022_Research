import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SOURCE = '~/Research/Spring2022_Research/ALMA_tight_binaries_Data_keepers2.csv'
df = pd.read_csv(SOURCE)
df = df.loc[:108] # non-data rows removed
J = df.loc[:,'J']
H = df.loc[:,'H']
K = df.loc[:,'K']
j_h = lambda h, k: 0.58 * (h-k) + 0.52
classII = (J-H) > j_h(H,K)
df2 = df[classII]
plt.plot(H-K, J-H, 'o')
plt.plot(H-K, j_h(H,K), '-')
ISR_vector = lambda x, c: 1.83*x + c
x1 = np.linspace(1.02/1.25,1.2,2)
x2 = np.linspace(0.618/1.25,1.2,2)
plt.plot(x1, ISR_vector(x1, -0.5),'-')
plt.plot(x2, ISR_vector(x2, -0.098),'-')
plt.show()
