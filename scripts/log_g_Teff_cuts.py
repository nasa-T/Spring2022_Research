import pandas as pd
import matplotlib.pyplot as plt

SOURCE = '~/Research/Spring2022_Research/ALMA_tight_binaries_Data.csv'
DESTINATION_KEEPERS = '~/Research/Spring2022_Research/ALMA_tight_binaries_Data_keepers.csv'
DESTINATION_ANNIHILATED = '~/Research/Spring2022_Research/ALMA_tight_binaries_Data_exterminated.csv'
df = pd.read_csv(SOURCE) # dataframe to manipulate
log_g = df.loc[:,'log(g) (cm/s2)']
teff = df.loc[:,'Teff']
cut_log_g = [1.7, 4.5, 5.0, 4.8] # log_g coords of points of cutoff
cut_teff = [3300, 4875, 3950, 2000] # Teff coords of points of cutoff
plt.plot(teff, log_g, 'o') # plot all data
plt.plot(cut_teff, cut_log_g, 'k-') # plot cutoff
plt.show()

m1 = 2.8/1575 # slope of line 1
b1 = 1.7 - m1 * 3300 # y-intercept of line 1
m2 = -1/1850 # slope of line 2
b2 = 4.5 - m2 * 4875 # y-intercept of line 2
m3 = 1/9750
b3 = 5.0 - m3 * 3950
y1 = lambda x1 : m1 * x1 + b1 # line between (1.7, 3300) and (4.5, 4875)
y2 = lambda x2 : m2 * x2 + b2 # line between (4.5, 4875) and (5.0, 3950)
y3 = lambda x3 : m3 * x3 + b3 # line between (5.0, 3950) and (4.8, 2000)

yep1 = df.loc[:,'log(g) (cm/s2)'] > y1(df.loc[:,'Teff'])
yep2 = df.loc[:,'log(g) (cm/s2)'] < y2(df.loc[:,'Teff'])
yep3 = df.loc[:,'log(g) (cm/s2)'] < y3(df.loc[:,'Teff'])
no1 = df.loc[:,'log(g) (cm/s2)'] < y1(df.loc[:,'Teff'])
no2 = df.loc[:,'log(g) (cm/s2)'] > y2(df.loc[:,'Teff'])
no3 = df.loc[:,'log(g) (cm/s2)'] > y3(df.loc[:,'Teff'])
stay1 = df.loc[:,'log(g) (cm/s2)'].isnull()
stay2 = df.loc[:,'Teff'].isnull()

df2 = df[yep1 & yep2 & yep3 | stay1 | stay2] # another dataframe for YSOs (made the cut)
df3 = df[-(yep1 & yep2 & yep3 | stay1 | stay2)] # objects that didn't make the cut
log_g = df2.loc[:,'log(g) (cm/s2)'] # coords for YSOs
teff = df2.loc[:,'Teff']
plt.plot(teff, log_g, 'o')
plt.plot(cut_teff, cut_log_g, 'k-')
plt.show()
print(df.loc[:,['Teff','log(g) (cm/s2)']])
print(df2.loc[:,['Teff','log(g) (cm/s2)']])
print(df3.loc[:,['Teff','log(g) (cm/s2)']]) 
df2.to_csv(DESTINATION_KEEPERS, index=False)
df3.to_csv(DESTINATION_ANNIHILATED, index=False)
# assume missing data is part of sample
# 3 sigma on parallax cuts
# check parallax-exterminated for log(g)-Teff

