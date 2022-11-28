from ggchemlib import *
import numpy as np
import matplotlib 
matplotlib.use("Agg") 
import pylab as pl
from scipy.interpolate import interp1d
from scipy.special import erfc

### The following ordered species are used for plotting:
order_TMC=['CO', 'H2O' ,'HCO+', 'H2CO','C3O', \
     'CH','C2H','C4H','C5H','C6H','C8H',\
     'NH3','N2H+','CN','C3N','C5N','HCN','HNC',\
     'HNC3','HC3N','HC5N','HC7N','HC9N',\
     'NO','CS','C2S','C3S','H2CS','HCS+',\
     'OH','OCS','CH3C5N',\
     'NS','CH4O','CH3C4H','O2','SO','SO2','H2S','HNCO']

def get_kappa(M, O, nH, ggchem=True):
    """
    calculate kappa according to the funtion defined by Garrod+2007.
    INPUT: 
    M      --- model file from ggchempy
    O      --- observations in dict, e.g. O[spec]=abun
    nH     --- gas density for the model
    ggchem --- True or False.
    
    OUTPUT: 
    kappa    --- kappa as function of time
    bestyear --- the best year
    """
    kappa=0.0
    N=0
    for key in O.keys():
        N+=1
        if ggchem==True: 
            X=M[key]/nH
        else:
            X=M[key]
        kappa = kappa + erfc(np.abs(np.log10(X)-np.log10(O[key]))/np.sqrt(2))
    kappa = kappa/float(N)
    maxkappa = np.max(kappa)
    ibestyear = np.argwhere(kappa==maxkappa)[0]
    if ggchem==True:
        bestyear = M['time'][ibestyear]
    else:
        bestyear = M['TIME'][ibestyear]
    return kappa, bestyear[0]

def compare(model, obs, nH, fn):    
    """
    Do the comparison between observations and models of TMC-1.
    INPUT:
    model ---- the model file of TMC1
    nH   ---- gas density used for the model
    fn   ---- output figure name
    """
    kappa, byr = get_kappa(model, obs, nH)
    time= model['time']
    sps = model.keys()
    
    pl.figure(figsize=(5,4))
    pl.subplots_adjust(bottom=0.15,right=0.95,top=0.95)
    pl.plot(time, kappa, 'k-')
    y12=pl.gca().get_ylim()
    pl.vlines(byr, y12[0],y12[1], linestyle='-')
    pl.xscale('log')
    pl.xlabel('Time (yr)')
    pl.ylabel('Mean confidence')
    pl.minorticks_on()
    pl.savefig(fn+'-kappa.pdf')
    pl.close('all')
    
    xtmc, xtmc0=[], []
    idx=[]
    i=0
    for sp in order_TMC:
        i+=1
        xtmc.append(obs[sp])
        xtmc0.append(interp1d(time, model[sp]/nH)(byr))
        idx.append(i)
        
    idx=np.array(idx)
    pl.figure(figsize=(7,12))
    pl.subplots_adjust(left=0.12,top=0.97,right=0.97,bottom=0.05)
    pl.barh(idx,     xtmc , height=0.7, align='center',color='none', edgecolor='k',label='Observations')
    pl.barh(idx+0.0, xtmc0, height=0.4, align='center',color='#A9A9A9', edgecolor='none', label=r'Models')
    pl.legend(loc=(0.5,0.5),ncol=1,frameon=False,fontsize=14)
    pl.xscale('log')
    y_order_TMC=[]
    for i in range(0,np.size(order_TMC),1):
        y_order_TMC.append( latex_species(order_TMC[i]) )
    pl.yticks(idx, y_order_TMC,rotation=0)
    #pl.xticks([],[])
    pl.xlim(1e-13,1.0e-2)
    pl.ylim(0,41)
    pl.tick_params(axis='y',which='both',top=False,bottom=False)
    pl.xlabel('Abundance',fontsize=15)
    pl.grid('on',axis='x')
    pl.savefig(fn+'.pdf')
    pl.close('all')
    
    print('The best year is %9.2e (yr)'%(byr))

if __name__=="__main__":
    ### load observations into a dict:
    ### Observed aundances from Agundez & Wakelam, 2013. 
    ### title: Chemistry of Dark Clouds: Databases, Networks, and Models 
    TMC={}
    spec_obs, abun_obs = np.loadtxt('data/TMC-1.obs', usecols=[0,1], skiprows=2, dtype='str', unpack=True)
    for i in range(spec_obs.size):
        TMC[spec_obs[i]]= abun_obs[i].astype(float)
    
    ### define the gas density (nH) for
    nH     = 2e4
    
    ### load a model:    
    model  = loadgg('../out/TMC1-updated.dat',skiprows=28)   
    compare(model, TMC, nH, 'tmc-1')
    
    ### load a model:    
    model_updated  = loadgg('../out/TMC1-updated.dat',skiprows=28)   
    compare(model_updated, TMC, nH, 'tmc-1-updated')
