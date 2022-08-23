#!/usr/bin/env python
from astropy.table import Table
import fitsio
import time
import astropy.io.fits as fits

__author__ = "Demitri Muna"
__date__   = "2017.05.09"

##  5.4G  survey-dr3-specObj-dr13.fits
# Download file from: http://portal.nersc.gov/project/cosmo/data/legacysurvey/dr3/external/survey-dr3-specObj-dr13.fits
filename = 'survey-dr3-specObj-dr13.fits'
columns = ['OBJID']

# Each test prints the first ten values of the column to guarantee
# the file has actually been read.

# File structure:
#
#  file: survey-dr3-specObj-dr13.fits
#  mode: READONLY
#  extnum hdutype         hduname[v]
#  0      IMAGE_HDU       
#  1      BINARY_TBL      DECALS-SDSS

def fitsio_read_1(filename):
	#
	# This goes straight to the table column, no messing about.
	#
	start_time = time.time()
	data = fitsio.read(filename, columns=['objid'], ext=1)
	print(data[0:10]) # print first ten values
	print("fitsio_read_1: {0} seconds".format(time.time() - start_time))

def fitsio_read_2(filename):
	#
	# This opens the file with more flexibility to read other data.
	#
	start_time = time.time()
	hdu_list = fitsio.FITS(filename) # open the file for reading
	decals_hdu = hdu_list[1] # select the HDU in question
	objid_col = decals_hdu.read(columns=['objid'])
	print(objid_col[0:10])
	print("fitsio_read_2: {0} seconds".format(time.time() - start_time))

def astropy_fits_read(filename):
	#
	# Setting mmap=True makes a big difference for large files.
	#
	start_time = time.time()
	hdu_list = fits.open(filename, memmap=True) # open the file for reading
	decals_hdu = hdu_list[1] # select the HDU in question - type:astropy.io.fits.hdu.table.BinTableHDU
	decals_table = decals_hdu.data
	print(decals_table['objid'][0:10])
	print("astropy_fits_read (mmap=True): {0} seconds".format(time.time() - start_time))

def astropy_table_read(filename):
	#
	# This works, but I don't like it as there is some "magic".
	# I never specified the HDU, and I don't like code making
	# assumptions like this. I don't know this class to be
	# able to make the code more clear/unambiguous.
	# It also reads the *full* table, so it will be slowest.
	#
	start_time = time.time()
	table = Table.read(filename, format="fits")
	print(table['OBJID'][0:10])
	print("astropy_table_read: {0} seconds".format(time.time() - start_time))

fitsio_read_1(filename)
fitsio_read_2(filename)
astropy_fits_read(filename)
astropy_table_read(filename)

# Times on my computer:
#   fitsio_read_1      - 2.2s
#   fitsio_read_2      - 2.2s
#   astropy_fits_read  - 4.0s
#   astropy_table_read - 36.2s
# 
# Take times to be approximate/representative - use timeit for more exhaustive testing.
#