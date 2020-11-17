from xspec import *
import os
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
from astropy.time import Time


def single_fit(x_pha, en=[0.5, 10.], plot=True):
    Fit.statMethod = "chi"
    Fit.nIterations = 100000
    Xset.parallel.error = 16
    Xset.parallel.leven = 16
    Xset.abund = 'wilm'
    Xset.xsect = 'vern'
    model = "tbabs*(diskbb+powerlaw+gauss)"
    pars = {1: 0.8, 2: 0.5, 3: 1000., 4: 2.0, 5: 0.5, 6: 6.5, 7: 1., 8: 0.1}
    fn = model.replace('*', '_')
    logfile = fn + '.log'
    Xset.openLog(logfile)
    s1 = Spectrum(x_pha)
    s1.ignore("0.-%.1f, %.1f-**"%(en[0],en[1]))
    AllData.ignore("bad")
    m1 = Model(model, setPars=pars)
    m1(6).values = [6.5, 0.01, 6.0, 6.0, 7., 7.]
    i1 = m1.nParameters
    i0 = m1.startParIndex
    compt = m1.componentNames
    Fit.perform()
    chi, dof = Fit.testStatistic, Fit.dof
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
    np.savez(zfile, model=model, phafile=x_pha, compt=compt, xx=xx, xx_err=xx_err, yy=yy, yy_err=yy_err, yy_fit=yy_fit, chisq=chi, dof=dof,
             par_name=par_name, par_value=par_value, par_error=par_error, par_range=par_range)
    AllData.clear()
    return zfile


def cal_flux(zfile, x_pha, Emin, Emax, en=[0.5, 10.], cal_error=True,flux_ini=[0.,0.]):
    tb=np.load(zfile)
    par_value=tb['par_value']
    model = "tbabs*(cflux*diskbb+cflux*powerlaw+gauss)"
    if flux_ini[0]==0:
        pars={1:par_value[0],2:Emin,3:Emax,4:-7.8,5:par_value[1], 6:par_value[2],7:Emin,8:Emax,9:-7.8,10:par_value[3],11:par_value[4],12:par_value[5],13:par_value[6],14:par_value[7]}
    else:
        pars={1:par_value[0],2:Emin,3:Emax,4:flux_ini[0],5:par_value[1], 6:par_value[2],7:Emin,8:Emax,9:flux_ini[1],10:par_value[3],11:par_value[4],12:par_value[5],13:par_value[6],14:par_value[7]}   
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
    s1.ignore("0.-%.1f, %.1f-**"%(en[0],en[1]))
    AllData.ignore("bad")

    m1 = Model(model,setPars=pars)
    m1(1).values = [2.5,0.01,1.5,1.5,4.0,4.0]
    #m1(5).values=[par_value[1],0.001,0.001,0.001,3.0,3.0]
    m1(12).values = [par_value[5], 0.01, 6.0, 6.0, 7.0, 7.0]
    m1(5).frozen,m1(10).frozen=True,True
    m1(6).frozen,m1(11).frozen=True,True
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
