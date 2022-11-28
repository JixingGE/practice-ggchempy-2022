def loadgg(ggfile,skiprows=24):
    import numpy as np
    f=open(ggfile,'r')
    d=f.readlines()
    x=d[skiprows].split()
    sp=x[:]
    f.close()
    
    datax=np.loadtxt(ggfile,unpack=True,skiprows=skiprows+1,dtype='str')
    ns = np.shape(datax)
    data = np.zeros(ns)
    for i in range(0,ns[0],):
        for j in range(0,ns[1],1):
            if 'e' not in datax[i][j] and 'E' not in datax[i][j]:
                data[i][j] = float(datax[i][j].replace('-','e-'))
            else:
                data[i][j] = float(datax[i][j])
    ab={}
    for i in range(0,np.size(sp),1):
        ab[sp[i]]=data[i][:]
    return ab
def latex_species(sp):
    xsp=r'${\rm '
    for x in sp:
        if x in '23456789':
            x='_'+x
        elif x in '+_':
            x='^{'+x+'}'
        xsp+=x
    xsp+='}$'
    if '[1_3C]' in xsp:
        ysp=xsp.replace('[1_3C]','^{13}C')
        return ysp
    else:
        return xsp
    
def latex_reaction(rn):
    xreaction=r''
    ir = 0
    irn=[1,1,0,1,1,0,0,0]
    for s in rn:
        ir+=1
        if s.strip()!='':
            if s.strip()=='E': s='e'
            sx = latex_species(s.strip())
            if ir==1:    xreaction += sx.strip()
            elif 1<ir<=3:xreaction += '+'+sx.strip()
            elif ir==4:  xreaction += r'$\rightarrow$ '+sx.strip()
            elif 4<ir<=8:xreaction += '+'+sx.strip()
    return xreaction
