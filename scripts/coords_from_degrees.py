from astropy import units as u
from astropy.coordinates import SkyCoord
import pandas as pd

# ORIGIN = "~/Research/Spring2022_Research/Kounkel2019_data.csv"
# DESTINATION = "~/Research/Spring2022_Research/Kounkel2019_data_convertedCoords.csv"

def extractData(filename):
    df = pd.read_csv(filename)
    return df

# convert coords in degrees to hms and dms
def fromDegrees(ra, dec):
    c = SkyCoord(ra=ra*u.degree,dec=dec*u.degree)
    return c.to_string('hmsdms')

origin = extractData(ORIGIN)
# apparently SkyCoord can't take more than 10000 entries at once, so we split the data
hmsdmsFirst = fromDegrees(origin.loc[0:9999,'RAJ2000'],origin.loc[0:9999,'DEJ2000'])
hmsdmsLast = fromDegrees(origin.loc[10000:,'RAJ2000'],origin.loc[10000:,'DEJ2000'])
print(origin)
origin['coords'] = hmsdmsFirst+hmsdmsLast
print(origin.loc[:,'coords'])
origin.to_csv(DESTINATION, index=False)
