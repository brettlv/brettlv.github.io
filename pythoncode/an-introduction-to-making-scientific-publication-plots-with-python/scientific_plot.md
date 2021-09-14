# [an-introduction-to-making-scientific-publication-plots](https://towardsdatascience.com/an-introduction-to-making-scientific-publication-plots-with-python-ea19dfa7f51e)

# Import required packages
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from pylab import cm
from matplotlib.pyplot import MultipleLocator
import matplotlib.dates as mdates
from matplotlib.transforms import Transform
from matplotlib.ticker import (AutoLocator, AutoMinorLocator)
from matplotlib.dates import DateFormatter
#import matplotlib.cm as cm
import matplotlib.colors as colors
from collections import OrderedDict


import numpy as np
import pandas as pd



import os
import sys
from time import time
from scipy.stats.stats import pearsonr
from scipy.stats import spearmanr#

from datetime import datetime, date, time, timezone
from datetime import datetime
from datetime import timedelta
from astropy.time import Time
from astropy.io import ascii

def datetime2mjd(x):
    mjd_ref=59000
    # calculate the difference between datenum of mjd_ref and mjd_ref
    # it should be a constant
    mjd_minus_mdates_num=mdates.date2num(convert_xaxis_time(mjd_ref))-mjd_ref
    
    x=mdates.date2num(x)# datenum of x
    y = x - mjd_minus_mdates_num
    return y # mjd of x

def mjd2datetime(x):
    mjd_ref=59000
    # calculate the difference between datenum of mjd_ref and mjd_ref
    # it should be a constant
    mjd_minus_mdates_num=mdates.date2num(convert_xaxis_time(mjd_ref))-mjd_ref
    y= x + mjd_minus_mdates_num# datenum of x
    y= mdates.num2date(y) #date of x
    return y



def datenums2mjd(x):
    #x=mdates.date2num(x)
    mjd_ref=59000
    mjd_minus_mdates_num=mdates.date2num(convert_xaxis_time(mjd_ref))-mjd_ref
    y = x - mjd_minus_mdates_num
    return y

def mjd2numsdate(x):
    mjd_ref=59000
    mjd_minus_mdates_num=mdates.date2num(convert_xaxis_time(mjd_ref))-mjd_ref
    y= x + mjd_minus_mdates_num
    #y= mdates.num2date(y)
    return y

def convert_xaxis_mjd(time):
    return Time(time).mjd

def convert_xaxis_time(mjd):
    return Time(mjd,format='mjd').to_datetime()


def date2yday(x):
    """
        x is in matplotlib datenums, so they are floats.
        """
    y = x - mdates.date2num(datetime(2018, 1, 1))
    return y

def yday2date(x):
    """
        return a matplotlib datenum (x is days since start of year of 2018)
        """
    y = x + mdates.date2num(datetime(2018, 1, 1))
    return y




def convert_partial_year(numbers):
    datetimes=[]
    for number in numbers:
        year = int(number)
        d = timedelta(days=(number - year)*(365 + is_leap(year)))
        day_one = datetime(year,1,1)
        date = d + day_one
        datetimes.append(date)
    return datetimes


def is_leap(year):
    if not year%4 and  year%100 or not year%400:
        return True
    return False


def convert_mjd(times):
    timesmjd=[]
    for i in times:
        timesmjd.append(Time(i).mjd)
    return timesmjd


def convert_date(times):
    timesdate=[]
    for i in times:
        timesdate.append(Time(i,format='mjd').datetime)
    return timesdate

def convert_date_single(time):
    timedate=Time(time,format='mjd').datetime
    return timedate







from astropy.table import Table
from astropy import units as u
from astropy import constants
from astropy.coordinates import SkyCoord
from astropy.visualization import hist
from astroML.datasets import fetch_imaging_sample, fetch_sdss_S82standards
from astroML.crossmatch import crossmatch_angular
from astropy import config as _config
from astroquery.irsa import Irsa
#from astroquery.ned import Ned

from collections import OrderedDict
from adjustText import adjust_text

import matplotlib.font_manager as fm
# Collect all the font names available to matplotlib
font_names = [f.name for f in fm.fontManager.ttflist]
#print(font_names)
#import matplotlib.font_manager as fm
# Rebuild the matplotlib font cache
#fm._rebuild()

%matplotlib inline
%config InlineBackend.figure_format='svg'



# Edit the font, font size, and axes width
mpl.rcParams['font.family'] = 'Avenir'
plt.rcParams['font.size'] = 18
plt.rcParams['axes.linewidth'] = 2

# Generate 2 colors from the 'tab10' colormap
colors = cm.get_cmap('tab10', 2)

def set_ax_tick(ax):
    ax.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', top='on',)
    ax.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in', top='on')
    ax.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', right='on')
    ax.yaxis.set_tick_params(which='minor', size=5, width=2, direction='in', right='on')

def set_ax_locator(ax,xma=1,xmi=0.2,yma=1,ymi=0.2):
    ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(xma))
    ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(xmi))
    ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(yma))
    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(ymi))

# drop index reset_index
def drop_index(data):
    data=data.reset_index(drop=True)
    return data

# drop columns
def drop_columns(data,column):
    data=data.drop(column,1)
    return data

##############################################################################################################################################################################

from astropy.cosmology import FlatLambdaCDM,Planck13,Planck15,z_at_value
from astropy import units as u
from astropy.cosmology import LambdaCDM
cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)




##############################################################################################################################################################################
# Use numpy.loadtxt to import our data
wavelength, samp_1_abs, samp_2_abs = np.loadtxt('Absorbance_Data.csv', unpack=True, delimiter=',', skiprows=1)
# Plot the two sample absorbances, using previously generated colors
ax.plot(wavelength, samp_1_abs, linewidth=2, color=colors(0), label='Sample 1')
ax.plot(wavelength, samp_2_abs, linewidth=2, color=colors(1), label='Sample 2')
# Set the axis limits
ax.set_xlim(370, 930)
ax.set_ylim(-0.2, 2.2)
# Edit the major and minor tick locations
ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(100))
ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(50))
ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.5))
ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.25))
# Add the x and y-axis labels
ax.set_xlabel('Wavelength (nm)', labelpad=10)
ax.set_ylabel('Absorbance (O.D.)', labelpad=10)
# Add the x-axis label with λ for wavelength
ax.set_xlabel(r'$\mathregular{\lambda}$ (nm)', labelpad=10)
# Create new axes object by cloning the y-axis of the first plot
ax2 = ax.twiny()
# Edit the tick parameters of the new x-axis
ax2.xaxis.set_tick_params(which='major', size=10, width=2, direction='in')
ax2.xaxis.set_tick_params(which='minor', size=7, width=2, direction='in')

# Function to convert energy (eV) to wavelength (nm)
def E_to_WL(E):
    return [1240/i for i in E]

# Add ticks manually to energy axis
ax2.xaxis.set_major_locator(mpl.ticker.FixedLocator(E_to_WL(np.linspace(1.5, 3.0, 4))))
ax2.xaxis.set_minor_locator(mpl.ticker.FixedLocator(E_to_WL(np.linspace(1.4, 3.2, 19))))
# Add tick labels manually to energy axis
ax2.set_xticklabels(['1.5', '2.0', '2.5', '3.0'])
# Add energy axis label
ax2.set_xlabel('Energy (eV)', labelpad=10)
# Set energy axis limits
ax2.set_xlim(370, 930)
# Add legend to plot
ax.legend(bbox_to_anchor=(1, 1), loc=1, frameon=False, fontsize=16)

# Save figure
plt.savefig('Final_Plot.png', dpi=400, transparent=False, bbox_inches='tight')
# Show figure
plt.show()


##############################################################################################################################################################################
# Use pandas to read csv or xlsx

data_csv=pd.read_csv(datapath,header=None,delimiter=',')
data_excel=pd.read_excel(datapath,sheet_name='',header=0)


# data_refill_nan
data_refill_nan.replace(to_replace=r'^\s*$',value=np.nan,regex=True,inplace=True)
data_refill_nan=data_refill_nan[data_refill_nan['Name'].notnull()]
data_refill_nan=data_refill_nan.reset_index(drop=True)


fig = plt.figure(figsize=(4,4))
fig.subplots_adjust(hspace=0.0, wspace = 0.0)
ax = fig.add_subplot(111)

ax.errorbars(data[x_label],data[y_label],xerr=,yerr=,
             color=color,
             marker=marker,
             label=label,
             )

texts = [ax.text(data[x_label][i], data[y_label][i]),
                 data[name_label][i],fontsize=9,)
                 for i in range(len(data))]

adjust_text(texts,ax=ax,
            #arrowprops=dict(arrowstyle='->', color='red',lw=0.5),
            expand_text=(1.15,1.15),
            expand_points=(1.15,1.15),
            expand_objects=(1.15, 1.15),
            expand_align=(1.15, 1.15),
            autoalign='xy',
            #only_move={'points':'x', 'text':'x'}
            ) #使用adjust_text



handles, labels = ax.get_legend_handles_labels()
hdl=handles
hdl = [h[0] for h in handles]# remove the errorbars
#use them in the legend
#by_label = OrderedDict(zip(labels, hdl)) #remove repeat label
labels_dict=dict(zip(labels, hdl)) #key,values
by_label=OrderedDict(sorted(labels_dict.items(),key=lambda t:t[0]))# sort based on labels
ax.legend(by_label.values(), by_label.keys(), bbox_to_anchor=(0.01, 0.99),
            loc=2, numpoints=1,ncol=1,fontsize=11.)




bottom,top=ax.set_ylim()
ax.set_ylim(top,bottom)
fig.savefig(imgpath,dpi=400, transparent=False, bbox_inches='tight')
##############################################################################################################################################################################




# Create figure object and store it in a variable called 'fig'
fig = plt.figure(figsize=(3, 3))
# Add axes object to our figure that takes up entire figure
ax = fig.add_axes([0, 0, 1, 1])
#Blank figure with axes at (0, 0) with a width and height of 1

# Add two axes objects to create a paneled figure
ax1 = fig.add_axes([0, 0, 1, 0.4])
ax2 = fig.add_axes([0, 0.6, 1, 0.4])
#Blank figure with two paneled axes, one at (0, 0) and the other at (0, 0.6), with widths of 1 and heights of 0.4

# Add two axes objects to create an inset figure
ax1 = fig.add_axes([0, 0, 1, 1])
ax2 = fig.add_axes([0.5, 0.5, 0.4, 0.4])
#Blank figure with inset — first axes takes up whole figure and second is at (0.5, 0.5) with a width and height of 0.4

#Remove Spines
# Hide the top and right spines of the axis
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

# Edit the major and minor ticks of the x and y axes
ax.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', top='on')
ax.xaxis.set_tick_params(which='minor', size=7, width=2, direction='in', top='on')
ax.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', right='on')
ax.yaxis.set_tick_params(which='minor', size=7, width=2, direction='in', right='on')
##############################################################################################################################################################################
