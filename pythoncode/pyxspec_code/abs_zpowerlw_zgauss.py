from xspec import *
import os
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
from astropy.time import Time


def single_fit(x_pha,obsid,filepath_root, en=[0.5, 8.], plot=True):
    Fit.statMethod = "chi"  #cstat
    Fit.nIterations = 100000
    Xset.parallel.error = 16
    Xset.parallel.leven = 16
    Xset.abund = 'wilm'
    Xset.xsect = 'vern'
    model = "tbabs*(zpowerlw+zgauss)"
    pars = {1: 0.0243, 2: 2.0, 3: 0.042436,4:1.0,
            5:  6.4, 6: 0.1, 7: 0.042436,   8:0.1}
    #pars = {1: 0.8, 2: 0.5, 3: 1000., 4: 2.0, 5: 0.5, 6: 6.5, 7: 1., 8: 0.1}
    # nh zgamma z zpnorm zgLineE zgSigma z zgnorm
    
    fn = model.replace('*', '_')
    logfile = fn + '.log'
    Xset.openLog(logfile)
    s1 = Spectrum(x_pha)

#corrarf=os.path.join(filepath_root,obsid,'%s_spec.corr.arf'%obsid)
#s1.response.arf=corrarf


    #s1.response.arf='src.arf'

    #s1.background='bkg.pha'
    #s1.response.arf='src.arf'
    #s1.response.rmf='src.rmf'
    s1.ignore("0.-%.1f, %.1f-**"%(en[0],en[1]))
    AllData.ignore("bad")
    m1 = Model(model, setPars=pars)
    m1(5).values = [6.5, 0.01, 6.0, 6.0, 7.0, 7.0]#Line
    #m1(5).values =6.4
    #m1(5).frozen=True
    m1(6).values = [0.1, 0.01, 1e-3, 1e-2, 0.001, 1.0] #Sigma
    #m1(6).values =0.1
    #m1(6).frozen=True

    m1(1).frozen=True



    i1 = m1.nParameters
    i0 = m1.startParIndex
    compt = m1.componentNames
    Fit.perform()
    chi, dof = Fit.testStatistic, Fit.dof
    chisq=chi/dof
    # if chi / dof < 2:
    #    Fit.error("%.f-%.f" % (i0, i1))
    #AllModels.calcFlux("%.1f %.1f err"%(en[0],en[1]))
    #flux, flux_lo, flux_up, ph, ph_lo, ph_up = s1.flux
    if plot:
        pfile = fn + '.ps'
        Plot.device = pfile + "/vcps"
        #Plot.device = "/xw"
        Plot.xAxis = "KeV"
        Plot.xLog = True
        Plot.yLog = True
        Plot.add = True
        #Plot.background = True
        Plot.setRebin(minSig=5., maxBins=20)
        Plot("data,delchi")
        #Plot("data resid")
        # Plot("ufspec")
        xx = Plot.x()
        xx_err = Plot.xErr()
        yy = Plot.y()
        yy_err = Plot.yErr()
        yy_fit = Plot.model()
        Plot.device = 'none'
    else:
        xx, yy, xx_err, yy_err, yy_fit = 0, 0, 0, 0, 0
    par_name = np.chararray(i1, itemsize=10)
    par_value = np.empty(i1)
    par_error = np.empty(i1)
    par_range = np.empty((2, i1))
    for k in range(i1):
        par_name[k] = m1(k + 1).name
        par_value[k] = m1(k + 1).values[0]
        par_range[:, k] = m1(k + 1).error[:-1]
        par_error[k] = m1(k + 1).sigma
    if os.path.isfile(fn + '.xcm'):
        os.remove(fn + '.xcm')
    Xset.save(fn, info='a')
    Xset.closeLog()
    zfile = fn + '.npz'
    print(par_name)
    np.savez(zfile, model=model, phafile=x_pha, compt=compt, xx=xx, xx_err=xx_err, yy=yy, yy_err=yy_err, yy_fit=yy_fit, chisq=chisq, dof=dof,
             par_name=par_name, par_value=par_value, par_error=par_error, par_range=par_range)
    AllData.clear()
    return zfile


def cal_flux(zfile,obsid,filepath_root, x_pha, Emin, Emax, en=[0.5, 8.0], cal_error=True,flux_ini=[0.]):
    tb=np.load(zfile)
    par_value=tb['par_value']
    model = "tbabs*(cflux*(zpowerlw+zgauss))"
    # 0,1,2,3,4,5,6,7 par_name
    #[b'nH' b'PhoIndex' b'Redshift' b'norm' b'LineE' b'Sigma' b'Redshift' b'norm']
    if flux_ini[0]==0:
        pars={1:par_value[0],2:Emin,3:Emax,4:-11.0,
        5:par_value[1], 6:par_value[2],7:par_value[3],8:par_value[4],
        9:par_value[5],10:par_value[6],11:par_value[7]}
    else:
        pars={1:par_value[0],2:Emin,3:Emax,4:flux_ini[0],
        5:par_value[1],6:par_value[2],7:par_value[3],8:par_value[4],
        9:par_value[5],10:par_value[6],11:par_value[7]}
    #nH,Emin,Emax,logF,
    #PhoIndex,Redshift,pnorm,LineE
    #Sigma,Redshift,norm

    Fit.statMethod = "chi"
    Fit.nIterations = 100000
    Xset.parallel.error = 16
    Xset.parallel.leven = 16
    Xset.abund = 'wilm'
    Xset.xsect = 'vern'

    fn = model.replace('*','_')
    logfile = fn+'.log'
    Xset.openLog(logfile)
    s1 = Spectrum(x_pha)
    #corrarf=os.path.join(filepath_root,obsid,'%s_spec.corr.arf'%obsid)
    #s1.response.arf=corrarf
    
    s1.ignore("0.-%.1f, %.1f-**"%(en[0],en[1]))
    AllData.ignore("bad")

    m1 = Model(model,setPars=pars)
    m1(1).frozen=True#nH
    #m1(1).values = [2.5,0.01,1.5,1.5,4.0,4.0]
    #m1(5).values=[par_value[1],0.001,0.001,0.001,3.0,3.0]
    m1(8).values = [par_value[4], 0.01, 6.0, 6.0, 7.0, 7.0] #LineE
    m1(9).values = [par_value[5], 0.01, 0.001, 0.001, 1.0, 1.0] #Sigma
    #m1(8).frozen=True
    m1(9).frozen=True,
   
    m1(5).values=[par_value[1],0.01,1.0,1.0,3.0,3.0]
    m1(5).frozen=True, #Gamma
    m1(6).frozen=True,
    
#m1(11).frozen=True,

    #AllModels.show()
    Fit.bayes='cons'
    Fit.perform()
    chi,dof = Fit.testStatistic, Fit.dof
    #if (chi/dof < 2) & cal_error:
    #    Fit.error("%.f-%.f"%(1,48))
    i1 = int(m1.nParameters)
    par_name=np.chararray(i1,itemsize=10)
    par_value=np.empty(i1)
    par_error = np.empty(i1)
    par_range=np.empty((2,i1))
    for k in range(i1):
        par_name[k] = m1(k+1).name
        par_value[k] = m1(k+1).values[0]
        par_range[:,k] = m1(k+1).error[:-1]
        par_error[k] = m1(k+1).sigma

    if os.path.isfile(fn+'.xcm'): os.remove(fn+'.xcm')
    Xset.save(fn, info='a')
    Xset.closeLog()
    zfile = fn+'.npz'
    np.savez(zfile,chisq=chi,dof=dof,par_name=par_name,par_value=par_value,par_error=par_error,par_range=par_range)
    AllData.clear()
    return zfile
