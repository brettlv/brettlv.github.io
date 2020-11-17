from xspec import *
import os
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
from astropy.time import Time

from abs_powerlaw import *
#from abs_zpowerlw_zgauss import *


def get_obsids(path):
    dirname=os.listdir(path)
    obsids=[]
    for i in dirname:
        if i.isdigit():
            obsids.append(i)
    obsids.sort()
    return obsids


def get_x_pha(filepath_root,obsid):
    return os.path.join(filepath_root,obsid),os.path.join(filepath_root,obsid,'%s_pileup_test_grp.pi'%obsid)


filepath_root='/Users/brettlv/CIAO/Mrk1018chandradata/chandra/'
obsids=get_obsids(filepath_root)

for i in range(1,len(obsids)):#len(obsids)
    obsid=obsids[i]
    obsid_pathdir,x_pha=get_x_pha(filepath_root,obsid)
    print(x_pha)
    os.chdir(obsid_pathdir)
    zfile=single_fit(x_pha,obsid,filepath_root, en=[0.5, 8.], plot=True)
    tb=np.load(zfile)
    dof=tb['dof']
    chisq=tb['chisq']
    gamma=tb['par_value'][1]
    gamma_error=tb['par_error'][1]
    model=tb['model']

print('current net rate/','net rate variance/','total rate/','predicted model rate')
for i in range(1,len(obsids)):
    obsid=obsids[i]
    obsid_pathdir,x_pha=get_x_pha(filepath_root,obsid)
    #print(x_pha)
    os.chdir(obsid_pathdir)
    model = "tbabs*powerlaw"
    fn = model.replace('*', '_')
    logfile = fn + '.log'
    zfile = fn + '.npz'
    tb=np.load(zfile)
    dof=tb['dof']
    chisq=tb['chisq']
    print(tb['rate'],)
