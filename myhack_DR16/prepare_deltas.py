'''                                                                                                                                                                                                                                           
Read the fits files of DR16 delta fields for Stripe 82                                                                                                                                                                                        
Calculate ra, dec, z and delta                                                                                                                                                                                                               
Store that data in a new file to be read by the VoidFinder algorithm.                                                                                                                                                                         
                                                                                                                                            
'''

print('Prepare the fits file for DR16S82 delta fields')

from astropy.table import Table
from astropy.io import fits
import os
import numpy as np

import matplotlib as mpl

import matplotlib.pyplot as plt

import matplotlib
#matplotlib.use("TkAgg")                                                                                                                                                                                                                      
from os import listdir
from os.path import isfile, join

in_directory='/global/homes/i/ineslie/myhack_DR16/Delta_LYA/Delta/'
os.chdir(in_directory)

#############################################################################                                                                                                                                                                 
#read_fits function reads the fits file, returns an astropy.Table with columns
#ra, dec, z, delta if the data is in the Stripe 82.                                                                                                                                                                                                                            
#############################################################################

def read_fits(namelist=('delta-89.fits.gz')):

    ra=list()
    dec=list()
    z=list()
    delta=list()
    calculated=Table()

    ra_min=-43
    ra_max=45
    dec_min=-1.25
    dec_max=1.25
    lambda_ref= 1215.668 #angstrom, reference wavelength. This is the Lyman-alpha spectral line for H.                                                                                                                                        
    i=0

    for filename in namelist:
        data = fits.open(filename)

        for hdu_num in range(1,len(data)):
            if data[hdu_num].header['RA']*(180/np.pi) > 180:
                data[hdu_num].header['RA']=data[hdu_num].header['RA']*(180/np.pi)-360
            else:
                data[hdu_num].header['RA']=data[hdu_num].header['RA']*(180/np.pi)
            if data[hdu_num].header['DEC']*(180/np.pi) > 180:
                data[hdu_num].header['DEC']=data[hdu_num].header['DEC']*(180/np.pi)-360
            else:
                data[hdu_num].header['DEC']=data[hdu_num].header['DEC']*(180/np.pi)

           # if ra_min <= data[hdu_num].header['RA'] <= ra_max and dec_min <= data[hdu_num].header['DEC'] <= dec_max:

            lambda_obs=10**(Table(data[hdu_num].data)['LOGLAM'])
                #lambda_rf=lambda_obs/((Table(data[0].data)['Z'][i]+1)
            data[hdu_num].header['RA']=data[hdu_num].header['RA']+90 #to make life easier for VF :)                                                                                                                                       
                #data[hdu_num].header['DEC']=data[hdu_num].header['DEC']+90 #not needed                                                                                                                                                       

            z_add=(lambda_obs-lambda_ref)/lambda_ref
            z.extend(z_add)
            delta.extend(Table(data[hdu_num].data)['DELTA'])
            ra.extend(data[hdu_num].header['RA']*np.ones(len(z_add)))
            dec.extend(data[hdu_num].header['DEC']*np.ones(len(z_add)))

    RA=Table.Column(ra, name='ra')
    DEC=Table.Column(dec, name='dec')
    Z=Table.Column(z, name='z')
    DELTA=Table.Column(delta, name='delta')

    calculated.add_column(RA)
    calculated.add_column(DEC)
    calculated.add_column(Z)
    calculated.add_column(DELTA)

    return calculated

onlyfiles = [f for f in listdir(in_directory) if isfile(join(in_directory, f))]

print(len(onlyfiles))

#onlyfiles.remove('Untitled.ipynb')
#onlyfiles.remove('prepare_deltas.py')
#onlyfiles.remove('prepare_deltas.py~')
onlyfiles.remove('deltafields_RAadded90.fits')
onlyfiles.remove('quasars.fits')
#onlyfiles.remove('.ipynb_checkpoints')
#onlyfiles.remove('#prepare_deltas.py#')

print(len(onlyfiles))

prepared=read_fits(onlyfiles)

#For vstack, I am worried about possible additional [ ]                                                                                                                                                                                       
#I am trying with my method instead.                                                                                                                                                                                                          


print('Necessary data calculated.')

prepared.write('deltafields_RAadded90.fits', format='fits', overwrite=True)

print('I have written the file. :)')
























