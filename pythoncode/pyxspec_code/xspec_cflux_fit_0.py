from xspec import *
import os
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
from astropy.time import Time


from abs_powerlaw import *
#from abs_zpowerlw import *
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
    return os.path.join(filepath_root,obsid),os.path.join(filepath_root,obsid,'%s_spec_grp.pi'%obsid)


filepath_root='/Users/brettlv/CIAO/Mrk1018chandradata/chandra/'

obsids=get_obsids(filepath_root)
obsid=obsids[0]
obsid_pathdir,x_pha=get_x_pha(filepath_root,obsid)
print(x_pha)
os.chdir(obsid_pathdir)
zfile=single_fit(x_pha,obsid,filepath_root, en=[4, 8.], plot=True)
from xspec import *
import os
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
from astropy.time import Time


from abs_powerlaw import *
#from abs_zpowerlw import *
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
    return os.path.join(filepath_root,obsid),os.path.join(filepath_root,obsid,'%s_spec_grp.pi'%obsid)


filepath_root='/Users/brettlv/CIAO/Mrk1018chandradata/chandra/'

obsids=get_obsids(filepath_root)
obsid=obsids[0]
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
print(obsid,'chisq:',chisq,'dof:',dof,model)



zfile_c=cal_flux(zfile,obsid,filepath_root, x_pha, Emin=2, Emax=10, en=[0.5, 8.], cal_error=True,flux_ini=[-11.0])

tb_c=np.load(zfile_c)
par_value=tb_c['par_value']
par_error=tb_c['par_error']

for i,j in zip(tb['par_name'],tb['par_value']):
    print(i,j)

print(par_value[4],par_error[4],10**par_value[3],10**(par_value[3]+par_error[3])-10**par_value[3])

#print(gamma,gamma_error,10**par_value[3],10**(par_value[3]+par_error[3])-10**par_value[3])
print('******')







zfile_c=cal_flux(zfile,obsid,filepath_root, x_pha, Emin=2, Emax=10, en=[0.5, 8.], cal_error=True,flux_ini=[-11.0])

tb_c=np.load(zfile_c)
par_value=tb_c['par_value']
par_error=tb_c['par_error']

for i,j in zip(tb['par_name'],tb['par_value']):
    print(i,j)

print(par_value[4],par_error[4],10**par_value[3],10**(par_value[3]+par_error[3])-10**par_value[3])

#print(gamma,gamma_error,10**par_value[3],10**(par_value[3]+par_error[3])-10**par_value[3])
print('******')




