import numpy as np
import os
from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt

# specify filenames
old_fname=os.environ['EBOSS_DR16']+'/qso/DR16Q/DR16Q_v4.fits'
new_fname='DR16Q_S82.fits'

# open original DR14Q file
hdulist=fits.open(old_fname)

# transform to table format, to modify later on
table_hdu = Table.read(hdulist[1])

# select objects in Stripe 82
mask = (table_hdu['DEC']>-1.26) & (table_hdu['DEC']<1.26) & ( (table_hdu['RA']<70.0) |  (table_hdu['RA']>300.0) )
print(np.sum(mask),'objects in Stripe 82')

# modify HDU with only quasars in Stripe 82
hdulist[1] = fits.table_to_hdu(table_hdu[mask])

# write new catalog to file
hdulist.writeto(new_fname,overwrite=True)

plt.plot(table_hdu['RA'],table_hdu['DEC'],'o',label='DR14Q')
plt.plot(table_hdu[mask]['RA'],table_hdu[mask]['DEC'],'o',label='in Stripe 82')
plt.xlabel('RA')
plt.ylabel('Dec')
plt.legend()
plt.show()

