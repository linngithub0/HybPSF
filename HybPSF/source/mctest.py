
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import emcee
from multiprocessing import Pool
np.random.seed(123)

def log_likelihood_SPCA(theta,star,PCs,weight):
    tmp=np.zeros(star.shape,dtype='float')
    for i in range(theta.shape[0]):
        tmp+=theta[i]*PCs[i,:]
    return -0.5 * np.sum((star-tmp)**2*weight)

def log_priori_SPCA(theta):
    Np=theta.shape[0]
    if Np==2:
        C0,C1=theta
        if -2.<C0<2. and -1.e-1<C1<1.e-1 :
            return 0.00
        else:
            return -np.inf
    if Np==3:
        C0,C1,C2=theta
        if -2.<C0<2. and -1.e-1<C1<1.e-1 and -1.e-2<C2<1.e-2:
            return 0.00
        else:
            return -np.inf
    if Np==4:
        C0,C1,C2,C3=theta
        if -2.<C0<2. and -1.e-1<C1<1.e-1 and -1.e-2<C2<1.e-2 and -1.e-2<C3<1.e-2:
            return 0.00
        else:
            return -np.inf
    if Np==5:
        C0,C1,C2,C3,C4=theta
        if (-2.<C0<2. and -1.e-1<C1<1.e-1 and -1.e-2<C2<1.e-2 and -1.e-2<C3<1.e-2 and 
            -1.e-2<C4<1.e-2):
            return 0.00
        else:
            return -np.inf
    if Np==6:
        C0,C1,C2,C3,C4,C5=theta
        if (-2.<C0<2. and -1.e-1<C1<1.e-1 and -1.e-2<C2<1.e-2 and -1.e-2<C3<1.e-2 and 
            -1.e-2<C4<1.e-2 and -1.e-2<C5<1.e-2):
            return 0.00
        else:
            return -np.inf

def log_probability_SPCA(theta,star,PCs,weight):
    lp = log_priori_SPCA(theta);
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood_SPCA(theta,star,PCs,weight)


def Coeff_deliv(star,PCs,weight):
    print("in Coeff_deliv")
    Nb=PCs.shape[1]
    pos=np.random.rand(PCs.shape[1]*10, PCs.shape[1])
    nwalkers, ndim = pos.shape
    sampler = emcee.EnsembleSampler(
        nwalkers, ndim, log_probability_SPCA, args=(star,PCs,weight)
    )
    sampler.run_mcmc(pos, 5000, progress=True);
    flat_samples = sampler.get_chain(discard=1000, thin=15, flat=True)
    tmp_coeff=[]
    for inc in range(Nb):
        #mcmcs = np.percentile(flat_samples[:, inc], [16, 50, 84])
        #tmp_coeff.append(list(mcmcs[1]))
        print(np.median(flat_samples[:, inc]))
    return(tmp_coeff)

def Coeff_mcmc(stars,PCs,weight):
    stars=np.array(stars);PCs=np.array(PCs);weight=np.array(weight)
    #write_mult_fits("/home/lnie/data/mcs.fits",stars,PCs,weight)
    Nstar=stars.shape[0]
    Mp=stars.shape[1]
    Nb=PCs.shape[1]
    #pool=Pool(processes=Nstar)
    print(stars.shape,PCs.shape,weight.shape)
    coeff=[]
    for ic in range(Nstar):
        print(ic)
        result=Coeff_deliv(stars[ic],PCs[ic],weight[ic])
        #result=pool.apply_async(Coeff_deliv,args=(stars[ic],PCs[ic],weight[ic]))
        coeff.append(result)
    #coeff=np.array(coeff)
    return(coeff)

hdu=fits.open("/Users/linn/Documents/code/JWST/fits/mcs.fits")
stars=hdu[0].data;PCs=hdu[1].data;weight=hdu[2].data
a=Coeff_mcmc(stars,PCs,weight)



