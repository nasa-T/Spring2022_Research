import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

SOURCE = '~/Research/Spring2022_Research/ALMA_tight_binaries_Data_keepers2.csv'

df = pd.read_csv(SOURCE)
df = df.loc[:108] # non-data rows removed

W1 = df.loc[:,'W1']
W2 = df.loc[:,'W2']
W3 = df.loc[:,'W3']
W4 = df.loc[:,'W4']
fig, ax = plt.subplots(2,2)
def barPlotting(dataset,wiseBand):
    """Takes dataset (Pandas dataframe) and wiseBand of interest (string) and returns quantities of each category"""
    detected = dataset[dataset.loc[:,'Disk Mass err']!=-10]
    non_detected = dataset[dataset.loc[:,'Disk Mass err']==-10]
    # Create different categories for combinations of WISE and ALMA data
    both = detected[detected.loc[:,wiseBand].notna()]
    WnoA = non_detected[non_detected.loc[:,wiseBand].notna()]
    AnoW = detected[detected.loc[:,wiseBand].isna()]
    noWnoA = non_detected[non_detected.loc[:,wiseBand].isna()]
    heights = np.array([len(both),len(WnoA),len(AnoW),len(noWnoA)])
    return heights

heights1 = barPlotting(df,'W1')
heights2 = barPlotting(df,'W2')
heights3 = barPlotting(df,'W3')
heights4 = barPlotting(df,'W4')
heights = np.array([[heights1, heights2], [heights3, heights4]])
for i in range(2):
    for j in range(2):
        for n in range(0,4):
            ax[i][j].text(n,heights[j][i][n],str(round(heights[j][i][n]/len(df),3)*100)+'%',ha='center')
ax[0][0].bar(['both','WISE no ALMA', 'ALMA no WISE', 'No Data'], heights1)
ax[0][0].set_title('ALMA-W1 correlation')
ax[1][0].bar(['both','WISE no ALMA', 'ALMA no WISE', 'No Data'], heights2)
ax[1][0].set_title('ALMA-W2 correlation')
ax[0][1].bar(['both','WISE no ALMA', 'ALMA no WISE', 'No Data'], heights3)
ax[0][1].set_title('ALMA-W3 correlation')
ax[1][1].bar(['both','WISE no ALMA', 'ALMA no WISE', 'No Data'], heights4)
ax[1][1].set_title('ALMA-W4 correlation')
plt.show()

