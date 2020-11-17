import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
#%matplotlib inline
#%config InlineBackend.figure_format='svg'
from pylab import cm
import matplotlib.font_manager as fm
# Collect all the font names available to matplotlib
font_names = [f.name for f in fm.fontManager.ttflist]
# Edit the font, font size, and axes width
mpl.rcParams['font.family'] = 'Avenir'
plt.rcParams['font.size'] = 18
plt.rcParams['axes.linewidth'] = 2


fig=plt.figure(8,6)
ax=fig.add_s    



def set_ax_axis_tick():
    ax.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', top='on')
    ax.xaxis.set_tick_params(which='minor', size=7, width=2, direction='in', top='on')
    ax.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', right='on')
    ax.yaxis.set_tick_params(which='minor', size=7, width=2, direction='in', right='on')


def plot_test():
    fig=plt.figure()
    a = np.arange(100)/10.0
    b= np.sin(a)
    plt.plot(a,b)
    plt.axhline(y=0.5)
    print(a[:3])
    print(b[:3])
