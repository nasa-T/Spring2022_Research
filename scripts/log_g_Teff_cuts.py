import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl


SOURCE = '~/Research/Spring2022_Research/ALMA_tight_binaries_Data2.csv'
# DESTINATION_KEEPERS = '~/Research/Spring2022_Research/ALMA_tight_binaries_Data_keepers.csv'
# DESTINATION_ANNIHILATED = '~/Research/Spring2022_Research/ALMA_tight_binaries_Data_exterminated.csv'
df = pd.read_csv(SOURCE) # dataframe to manipulate
tooClose = df[df.loc[:,'p'] > 3.5]
tooFar = df[df.loc[:,'p'] < 2]
Goldilocks = df[(df.loc[:,'p'] >= 2) & (df.loc[:,'p'] <= 3.5)]
Glog_g = Goldilocks.loc[:,'log(g) (cm/s2)']
Gteff = Goldilocks.loc[:,'Teff']
Clog_g = tooClose.loc[:,'log(g) (cm/s2)']
Cteff = tooClose.loc[:,'Teff']
Flog_g = tooFar.loc[:,'log(g) (cm/s2)']
Fteff = tooFar.loc[:,'Teff']
cut_log_g = [1.7, 4.5, 5.0, 4.8] # log_g coords of points of cutoff
cut_teff = [3300, 4875, 3950, 2000] # Teff coords of points of cutoff
# print(df.loc[:,'p'])
norm = mpl.colors.Normalize(vmin=0.5, vmax=5)
cmap = mpl.cm.coolwarm
# plt.scatter(teff, log_g, c=df.loc[:,'p'], s=10, alpha=1, cmap=cmap)
# plt.plot(Fteff, Flog_g, 'ro')
# plt.plot(Gteff, Glog_g, 'bo')
# plt.plot(Cteff, Clog_g, 'go')
plt.errorbar(Fteff, Flog_g, fmt='ro', xerr=tooFar.loc[:,'Teff_error'], yerr=tooFar.loc[:,'log(g)_error'],markersize=3)
plt.errorbar(Gteff, Glog_g, fmt='bo', xerr=Goldilocks.loc[:,'Teff_error'], yerr=Goldilocks.loc[:,'log(g)_error'],markersize=3)
plt.errorbar(Cteff, Clog_g, fmt='go', xerr=tooClose.loc[:,'Teff_error'], yerr=tooClose.loc[:,'log(g)_error'],markersize=3)
plt.legend(['Too Far', 'Just Right', 'Too Close'])
plt.xlabel('teff')
plt.ylabel('log(g)')
# plt.errorbar(df.loc[:,'Teff'], df.loc[:,'log(g) (cm/s2)'], fmt='o', xerr=df.loc[:,'Teff_error'], yerr=df.loc[:,'log(g)_error'])
# plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), label='parallax', orientation='vertical')
# plt.plot(teff, log_g, 'o') # plot all data
plt.plot(cut_teff, cut_log_g, 'k-') # plot cutoff
plt.show()

m1 = 2.8/1575 # slope of line 1
b1 = 1.7 - m1 * 3300 # y-intercept of line 1
m2 = -1/1850 # slope of line 2
b2 = 4.5 - m2 * 4875 # y-intercept of line 2
m3 = 1/9750 # slope of line 3
b3 = 5.0 - m3 * 3950 # y-intercept of line 3
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

tooClose2 = df2[df2.loc[:,'p'] > 3.5]
tooFar2 = df2[df2.loc[:,'p'] < 2]
Goldilocks2 = df2[(df2.loc[:,'p'] >= 2) & (df2.loc[:,'p'] <= 3.5)]
Glog_g2 = Goldilocks2.loc[:,'log(g) (cm/s2)']
Gteff2 = Goldilocks2.loc[:,'Teff']
Clog_g2 = tooClose2.loc[:,'log(g) (cm/s2)']
Cteff2 = tooClose2.loc[:,'Teff']
Flog_g2 = tooFar2.loc[:,'log(g) (cm/s2)']
Fteff2 = tooFar2.loc[:,'Teff']

plt.plot(Fteff2, Flog_g2, 'ro')
plt.plot(Gteff2, Glog_g2, 'bo')
plt.plot(Cteff2, Clog_g2, 'go')
plt.legend(['Too Far', 'Just Right', 'Too Close'])

# plt.plot(teff, log_g, 'o')
plt.plot(cut_teff, cut_log_g, 'k-')
# plt.scatter(teff, log_g, c=df2.loc[:,'p'], s=10, alpha=1, cmap=cmap)
# plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), label='parallax', orientation='vertical')
plt.show()
# print(df.loc[:,['Teff','log(g) (cm/s2)']])
# print(df2.loc[:,['Teff','log(g) (cm/s2)']])
# print(df3.loc[:,['Teff','log(g) (cm/s2)']]) 
# df2.to_csv(DESTINATION_KEEPERS, index=False)
# df3.to_csv(DESTINATION_ANNIHILATED, index=False)
# assume missing data is part of sample
# 3 sigma on parallax cuts
# check parallax-exterminated for log(g)-Teff
# color bar for parallaxes

