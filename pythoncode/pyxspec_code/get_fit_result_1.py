from xspec import *
import os
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
from astropy.time import Time



def get_obsids(path):
    dirname=os.listdir(path)
    obsids=[]
    for i in dirname:
        if i.isdigit():
            obsids.append(i)
    obsids.sort()
    return obsids


def get_x_pha(filepath_root,obsid):
    return os.path.join(filepath_root,obsid),os.path.join(filepath_root,obsid,'%s_spec_grp.pi'%obsid)


filepath_root='/Users/brettlv/CIAO/Mrk1018chandradata/chandra/'
obsids=get_obsids(filepath_root)




for i in range(1,len(obsids)):#len(obsids)
    obsid=obsids[i]
    obsid_pathdir,x_pha=get_x_pha(filepath_root,obsid)
    #print(x_pha)
    os.chdir(obsid_pathdir)
    
    #model = "tbabs*(zpowerlw+zgauss)"
    model="tbabs*zpowerlw"
    fn = model.replace('*', '_')
    zfile = fn + '.npz'
    
    tb=np.load(zfile)
    dof=tb['dof']
    chisq=tb['chisq']
    gamma=tb['par_value'][1]
    gamma_error=tb['par_error'][1]


#model = "tbabs*(cflux*(zpowerlw+zgauss))"
    model="tbabs*(cflux*zpowerlw)"
    
    fn = model.replace('*', '_')
    zfile = fn + '.npz'
    tb=np.load(zfile)
    par_name=tb['par_name']
    #print(par_name)
    par_value=tb['par_value']
    par_error=tb['par_error']
    print(obsid,'chisq:','{:.2f}'.format(chisq),'dof:',dof)
    #print(gamma,gamma_error,10**par_value[3],10**(par_value[3]+par_error[3])-10**par_value[3])
    print('G:','{:.2f}'.format(gamma),'{:.2f}'.format(gamma_error),'F:','{:.2e}'.format(10**par_value[3]),'{:.2e}'.format(10**(par_value[3]+par_error[3])-10**par_value[3]))
    print()
#print('******')



