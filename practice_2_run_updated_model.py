from ggchempylib import ggchempy, loadgg, latex_species, run_GUI
import os
import numpy as np
gg = ggchempy

### The following three switches (logical) are only used for this script:
useGUI=False       ## If True, PyQt5 is needed.
### Main code:
if useGUI:
    run_GUI()
else:
    
    #################################### the static models:
    GGCHEM={}
    modelnames = ['TMC1-updated']

    ### define common parameters:
    gg.ggpars.d2gmr = 0.01
    gg.ggpars.ti = 1.0
    gg.ggpars.tf = 1.0e+9
    gg.ggpars.ggfiles= ['in/network2.txt','in/ed.txt','in/iabun.txt']
    gg.gas.Zeta = 1.3e-17
    gg.dust.surface.Rdb=0.77
    
    ### the five static models:
    pars={}        # nH,     Tgas, Av,   chi,  Tdust
    pars['TMC1-updated']   =[2.00e4, 10.0, 10.0, 1.0,  10.0]
    
    for model in modelnames:
        gg.gas.nH, gg.gas.T, gg.gas.Av, gg.gas.Chi, gg.dust.T = pars[model]
        GGCHEM[model] = gg.run(model)  ## diectory "out" will be created to save model. e.g. "out/TMC1.dat"
