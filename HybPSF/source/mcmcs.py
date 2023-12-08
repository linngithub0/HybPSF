import numpy as np
import numpy
import ctypes
import os
from astropy.io import fits
import emcee
import webbpsf
from scipy import interpolate
from multiprocessing import Pool

#from psf_fit import web_psf_rec,wgtmap,model2pos,size
#from C_tools_lib import centriod_psf

detector=dict({"nrca1":"nrca1","nrca2":"nrca2","nrca3":"nrca3","nrca4":"nrca4",
               "nrcb1":"nrcb1","nrcb2":"nrcb2","nrcb3":"nrcb3","nrcb4":"nrcb4",
               "nrcalong":"nrca5","nrcblong":"nrcb5"})


def gaus_estimate(image):
    nx=image.shape[0];ny=image.shape[1];#print("x=%f,y=%f"%(nx,ny))
    val=np.empty((0,1),float);
    r2min=(min(nx,ny)/2.)**2;#print("r2min=%f"%(r2min))
    for i in range(nx):
        for j in range(ny):
            x=i-nx/2.+0.5;
            y=j-ny/2.+0.5
            c=image[i][j]
            r2=x*x+y*y
            if r2>r2min and c!=0:
                val=np.append(val,[[c],],axis=0)
    return(np.mean(val),np.std(val))


def wgtmap(image,gain=1):
    Ng1=image.shape[0]
    Ng2=image.shape[1]
    weight=np.zeros(image.shape,dtype='float')
    mean,sigma=gaus_estimate(image)
    #print("method=%f,sigma=%f"%(mean,sigma))
    for i in range(Ng1):
        for j in range(Ng2):
            detx=sigma*sigma*gain*gain
            dety=gain*np.fabs(image[i][j])
            if dety<detx : error=sigma
            else : error=np.sqrt(dety+detx)/gain
            #weight[i][j]=1./(error**2)
            weight[i][j]=error
    return(weight)

def S2N(image,gain=1):
    Ng1=image.shape[0]
    Ng2=image.shape[1]
    weight=np.zeros(image.shape,dtype='float')
    mean,sigma=gaus_estimate(image)
    #print("method=%f,sigma=%f"%(mean,sigma))
    for i in range(Ng1):
        for j in range(Ng2):
            detx=sigma*sigma*gain*gain
            dety=gain*np.fabs(image[i][j])
            if dety<detx : error=sigma
            else : error=np.sqrt(dety+detx)/gain
            #weight[i][j]=1./(error**2)
            weight[i][j]=error
    snr=np.sum(image)/np.sqrt(np.sum(weight**2))
    return(snr)

def write_fits(fitsname,data,head1=[],head2=[]):
    wdata=data;
    if head1==None :
        #print("head1 none",fitsname)
        hdu = fits.PrimaryHDU(wdata)
        hdul = fits.HDUList([hdu])
        hdul.writeto(fitsname,overwrite=True)
    if type(head1) is np.ndarray and type(head2) is list:
        #print("head2 none",fitsname)
        hdu = fits.PrimaryHDU(wdata,header=head1)
        hdul = fits.HDUList([hdu])
        hdul.writeto(fitsname,overwrite=True)
    if type(head2) is np.ndarray:
        #print("head1 not none",fitsname)
        hdu = fits.PrimaryHDU(header=head1)
        hdu1= fits.ImageHDU(wdata,header=head2)
        hdul = fits.HDUList([hdu,hdu1])
        hdul.writeto(fitsname,overwrite=True)
        
def write_mult_fits(fitsname,data1=[],data2=[],data3=[],data4=[]):
    if type(data2) is list :
        hdu = fits.PrimaryHDU(data1)
        hdul = fits.HDUList([hdu])
        hdul.writeto(fitsname,overwrite=True) 
    if type(data2) is np.ndarray and type(data2) is list:
        hdu = fits.PrimaryHDU(data1)
        hdu1=fits.ImageHDU(data2)
        hdul = fits.HDUList([hdu,hdu1])
        hdul.writeto(fitsname,overwrite=True) 
    if type(data3) is np.ndarray and type(data4) is list:
        hdu = fits.PrimaryHDU(data1)
        hdu1=fits.ImageHDU(data2)
        hdu2=fits.ImageHDU(data3)
        hdul = fits.HDUList([hdu,hdu1,hdu2])
        hdul.writeto(fitsname,overwrite=True)
    if type(data4) is np.ndarray:
        hdu = fits.PrimaryHDU(data1)
        hdu1=fits.ImageHDU(data2)
        hdu2=fits.ImageHDU(data3)
        hdu3=fits.ImageHDU(data4)
        hdul = fits.HDUList([hdu,hdu1,hdu2,hdu3])
        hdul.writeto(fitsname,overwrite=True)
        
def size(image,center,sigma):
    '''
    size(image,center,sigma)
    center[0]=cx;center[1]=cy
    '''
    nx=image.shape[0];ny=image.shape[1]
    #print("nx=%d,ny=%d"%(nx,ny))
    W=0;R11=0;R22=0;R12=0;R2=0.;k=0;
    nh=(np.min([nx,ny])*0.5)**2
    #print("nh=%d,cx=%f,cy=%f"%(nh,center[0],center[1]))
    scale=0.5/(sigma**2)
    for i in range(nx):
        for j in range(ny):
            x=i-center[0];y=j-center[1]
            r2=x**2+y**2;weight=np.exp(-r2*scale)
            if r2<nh:
                W=W+image[i][j]*weight
                R11=R11+x*x*image[i][j]*weight
                R22=R22+y*y*image[i][j]*weight
                R12=R12+x*y*image[i][j]*weight
                #if k<10 :print(image[i][j],weight,sigma)
                #k=k+1;
    R11=R11/W;R22=R22/W;R12=R12/W
    #print("message:w11=%f,w12=%f,w22=%f,w=%f"%(R11,R12,R22,W))
    e1=(R11-R22)/(R11+R22)
    e2=(2.*R12)/(R11+R22)
    #print("message")
    R2=R11+R22
    #print(e1,e2,R2)
    if R2<=0 : R2=0.001
    if R2<=0 or np.fabs(e1)>=1. or np.fabs(e2)>=1.:
        #raise Exception("R2<=0\n")
        print("R2<=0\n");
        return(e1,e2,R2)
    else:
        return(e1,e2,R2)


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
        if -2.<C0<2. and -1.<C1<1. and -1.<C2<1.:
            return 0.00
        else:
            return -np.inf
    if Np==4:
        C0,C1,C2,C3=theta
        if (-2.<C0<2. and -1.<C1<1. and -1.<C2<1. and -1.<C3<1.):
            return 0.00
        else:
            return -np.inf
    if Np==5:
        C0,C1,C2,C3,C4=theta
        if ( -2.<C0<2. and -1.<C1<1. and -1.<C2<1. and -1.<C3<1. and 
            -1.<C4<1.):
            return 0.00
        else:
            return -np.inf
    if Np==6:
        C0,C1,C2,C3,C4,C5=theta
        if ( -2.<C0<2. and -1.<C1<1. and -1.<C2<1. and -1.<C3<1. and 
            -1.<C4<1. and -1.<C5<1.):
            return 0.00
        else:
            return -np.inf

def log_probability_SPCA(theta,star,PCs,weight):
    lp = log_priori_SPCA(theta);
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood_SPCA(theta,star,PCs,weight)


def Coeff_deliv(pos,star,PCs,weight):
    nwalkers, ndim = pos.shape
    sampler = emcee.EnsembleSampler(
        nwalkers, ndim, log_probability_SPCA, args=(star,PCs,weight)
        )
    sampler.run_mcmc(pos, 5000, progress=False);
    flat_samples = sampler.get_chain(discard=500, thin=1, flat=True)
    tmp_coeff=[]
    for inc in range(ndim):
        mcmcs = np.percentile(flat_samples[:, inc], [16, 50, 84])
        q = np.diff(mcmcs)
        tmp_coeff.append(mcmcs[1])
    return(tmp_coeff)

def Coeff_mcmc(stars,PCs,weight,coeffs):
    stars=np.array(stars);PCs=np.array(PCs);weight=np.array(weight);coeffs=np.array(coeffs)
    #write_mult_fits("mcs.fits",stars,PCs,weight,coeffs)
    Nstar=stars.shape[0]
    Mp=stars.shape[1]
    Nb=PCs.shape[1]
    pool=Pool(processes=Nstar)
    print(stars.shape,PCs.shape,weight.shape)
    coeff=[]
    results=[]
    for ic in range(Nstar):
        pos=np.random.randn(PCs.shape[1]*10, PCs.shape[1])*10**(-2);
        result=pool.apply_async(func=Coeff_deliv,args=(pos,stars[ic],PCs[ic],weight[ic],))
        results.append(result)
    for re in results:
        coeff.append(re.get())
    pool.close() 
    pool.join()
    return(coeff)


'''
def Coeff_mcmc(stars,PCs,weight,coeffs):
    stars=np.array(stars);PCs=np.array(PCs);weight=np.array(weight)
    write_mult_fits("/home/lnie/data/mcs.fits",stars,PCs,weight,coeffs)
    Nstar=stars.shape[0]
    Mp=stars.shape[1]
    Nb=PCs.shape[1]
    print(stars.shape,PCs.shape,weight.shape)
    coeff=[]
    for ic in range(Nstar):
        print("coeff ic",ic)
        pos=np.random.randn(PCs.shape[1]*10, PCs.shape[1])*10**(-2)
        print("pos",pos.shape)
        nwalkers, ndim = pos.shape
        sampler = emcee.EnsembleSampler(
            nwalkers, ndim, log_probability_SPCA, args=(stars[ic],PCs[ic],weight[ic])
            )
        sampler.run_mcmc(pos, 5000, progress=True);
        flat_samples = sampler.get_chain(discard=500, thin=15, flat=True)
        tmp_coeff=[]
        for inc in range(Nb):
            mcmcs = np.percentile(flat_samples[:, inc], [16, 50, 84])
            q = np.diff(mcmcs)
            tmp_coeff.append(mcmcs[1])
        coeff.append(tmp_coeff)
    #coeff=np.array(coeff);
    #print(coeff.shape)
    return(coeff)
'''




def File_judge(stamp_name):
    print("File_judge")
    #print(stamp_name)
    try:
        #print("0");os.system("ls "+stamp_name)
        hdu=fits.open(stamp_name);
        stamps=hdu[0].data;
        Nobj=stamps.shape[0];
        print(stamps.shape);
        if Nobj>6:
            return(1)
        else:
            if Nobj<2:return(-1)
            else:return(0)
    except:
        return(-1)


def Get_data(stamp_name,ofname,dirs,PCdirs,det,filte,Nim):
    print("use only webpsf")
    dirs=dirs+"NIRCam/"+filte+"/"
    PCdirs=PCdirs+"NIRCam/"+filte+"/"
    print(stamp_name,ofname,dirs,PCdirs,det,filte)
    tmpname=stamp_name.split("/");l0=len(tmpname);l1=len(tmpname[l0-1]);
    slicename=tmpname[l0-1][0:l1-len("_star.fits")]
    hdu=fits.open(ofname);OPD=hdu[0].header['DATE-BEG'];print("OPD:",OPD);
    PHOTUJA2=hdu['SCI'].header['PHOTUJA2'];XPOSURE=hdu['SCI'].header['XPOSURE']
    PHOTMJSR=hdu['SCI'].header['PHOTMJSR']
    PIXAR_A2=hdu['SCI'].header['PIXAR_A2']
    cat_name=dirs+det+'/'+"catalogue/star/"+slicename+"_star.cat"
    mask_name=dirs+det+'/'+"catalogue/star_stamps/"+slicename+"_mask.fits"
    cat=np.loadtxt(cat_name);print(cat)
    hdu=fits.open(stamp_name);instar=hdu[0].data*XPOSURE*PHOTMJSR;
    Nobj=instar.shape[0]  
    hdu=fits.open(mask_name)
    inmask=hdu[0].data;print("inmask shape:",inmask.shape)
    Ng0=instar.shape[1];Nh=int(Ng0/2);osam=2
    for ic in range(Nobj):
        instar[ic]=instar[ic].T
        inmask[ic]=inmask[ic].T
    spos=np.zeros((Nobj,2),dtype='float')
    wgt=np.zeros(instar.shape,dtype='float')
    for ic in range(Nobj):
        spos[ic][0]=cat[ic][1];spos[ic][1]=cat[ic][2]
        rx=cat[ic][1]+0.5;ry=cat[ic][2]+0.5
        xl=int(rx-Nh);xh=int(rx+Nh);yl=int(ry-Nh);yh=int(ry+Nh)
        wgt[ic]=inmask[ic][:,:]
    print("call webbpsf")
    imodel=np.zeros((Nobj,(Ng0+Nim)*osam,(Ng0+Nim)*osam),dtype='float')
    nrc = webbpsf.NIRCam()
    nrc.filter =filte;print(filte)
    nrc.detector=detector[det].upper();print(detector[det].upper())
    nrc.load_wss_opd_by_date(OPD);print(OPD)
    for ic in range(Nobj):
        print(slicename+" mixpsf imodel",ic+1,spos[ic][0],spos[ic][1])
        nrc.detector_position=(spos[ic][0],spos[ic][1])
        psf = nrc.calc_psf(fov_pixels=Ng0+Nim,oversample=2)
        imodel[ic]=(psf[2].data/np.sum(psf[2].data)).T
        #imodel[ic]/=np.sqrt(np.sum(imodel[ic])**2)
    write_fits("fits/cmodel.fits",imodel)
    instar=list(instar);spos=list(spos);wgt=list(wgt);imodel=list(imodel)
    for ic in range(Nobj):
        instar[ic]=list(instar[ic])
        spos[ic]=list(spos[ic])
        wgt[ic]=list(wgt[ic])
        imodel[ic]=list(imodel[ic])
        for j in range(Ng0):
            instar[ic][j]=list(instar[ic][j])
            wgt[ic][j]=list(wgt[ic][j])
        for j in range((Ng0+Nim)*osam):
            imodel[ic][j]=list(imodel[ic][j])
    return(instar,spos,wgt,imodel)


def log_likelihood_poly(theta,posxy,coeff,error):
    Ci=np.zeros(coeff.shape[0])
    for ic in range(coeff.shape[0]):
        Ci[ic]=np.sum(posxy[:,ic]*theta)
    return -0.5*np.sum((coeff-Ci)**2./error**2)


def log_priori_poly(theta):
    count=1
    for thetai in theta:
        if thetai>-0.1 and thetai<0.1:
            count*=1
        else:
            count*=0
    if count==1:
        return(0.0)
    else:
        return(-np.inf)

def log_probability_ploy(theta,posxy,coeff,error):
    lp = log_priori_poly(theta);
    if not np.isfinite(lp):
        return(-np.inf)
    chi2=log_likelihood_poly(theta,posxy,coeff,error)
    return(lp + chi2)

def mcmc_poly(pos,val,error,ipos,polyorder):
    Nstar=pos.shape[0];
    l0=0
    for i in range(polyorder+1):
        l0+=i+1
    #print("l0=",l0)
    xy=np.zeros((l0,Nstar),dtype='float')
    ival=np.zeros(ipos.shape[0],dtype='float')
    ixy=np.zeros((l0,ipos.shape[0]),dtype='float')
    poly_coeff=1e-2 * np.random.randn(l0*10, l0)
    for ic in range(Nstar):
        x=pos[ic][0];y=pos[ic][1]
        l=0
        for i in range(polyorder+1):
            for j in range(polyorder+1):
                if (i+j)<=polyorder:
                    xy[l][ic]=x**i*y**j 
                    l+=1
    for ic in range(ipos.shape[0]):
        x=ipos[ic][0];y=ipos[ic][1]
        l=0
        for i in range(polyorder+1):
            for j in range(polyorder+1):
                if (i+j)<=polyorder:
                    ixy[l][ic]=x**i*y**j 
                    l+=1
    nwalkers, ndim = poly_coeff.shape
    sampler = emcee.EnsembleSampler(
        nwalkers, ndim, log_probability_ploy, args=(xy,val,error))
    sampler.run_mcmc(poly_coeff, 5000, progress=True);
    flat_samples = sampler.get_chain(discard=100, thin=15, flat=True)
    theta=np.zeros(l0,dtype='float')
    for i in range(ndim):
        mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
        q = np.diff(mcmc)
        theta[i]=mcmc[1]
    for ic in range(ipos.shape[0]):
        ival[ic]=np.sum(theta*ixy[:,ic])
    return(ival)


def polyfit2d(pos,val,error,tar_pos,polyorder):
    """
    pos,star position
    val, value at star position
    error, errors of val
    tar_pos, target position
    polyorder
    """
    Nstar=pos.shape[0];iNstar=tar_pos.shape[0]
    l0=0
    for i in range(polyorder+1):
        l0+=i+1
    #print("l0=",l0)
    xy=np.zeros((l0,Nstar),dtype=np.float)
    a=np.zeros((l0,l0),dtype=np.float)
    b=np.zeros(l0,dtype=np.float)
    ival=np.zeros(iNstar,dtype=np.float)
    ixy=np.zeros((l0,iNstar),dtype=np.float)
    for ic in range(Nstar):
        x=pos[ic][0];y=pos[ic][1]
        l=0
        for i in range(polyorder+1):
            for j in range(polyorder+1):
                if (i+j)<=polyorder:
                    xy[l][ic]=x**i*y**j 
                    l+=1
    for l in range(l0):
        for f in range(l0):
            for ic in range(Nstar):
                a[l][f]+=xy[l][ic]*xy[f][ic]*(1./error[ic]**2)
        for ic in range(Nstar):
            b[l]+=xy[l][ic]*val[ic]*(1./error[ic]**2)
    poly=np.linalg.solve(a,b);#print(poly)
    for ic in range(iNstar):
        x=tar_pos[ic][0];y=tar_pos[ic][1]
        l=0
        for i in range(polyorder+1):
            for j in range(polyorder+1):
                if (i+j)<=polyorder:
                    ixy[l][ic]=x**i*y**j 
                    l+=1
    for ic in range(iNstar):
        ival[ic]=np.sum(poly*ixy[:,ic])
    return(ival)


def coeff2psf(spos,coeffs,PCs,gpos,degrees):
    """
    spos, star position
    coeffs, PC coeffs of stars
    PCs, PCs
    gpos, galaxy position
    degrees,
    """
    iNstar=gpos.shape[0];Nstar=spos.shape[0]
    Ng=PCs.shape[1]
    Nb=PCs.shape[0]
    icoeff=np.zeros((iNstar,Nb),dtype='float')
    for ic in range(Nb):
        #icoeff[ic]=mcmc_poly(spos,coeffs[ic,:,0],coeffs[ic,:,1],gpos,degrees)
        icoeff[:,ic]=polyfit2d(spos,coeffs[:,ic,0],coeffs[:,ic,1],gpos,degrees);#print(icoeff[:,ic])
    rPSF=np.zeros((iNstar,Ng,Ng),dtype='float')
    for ic in range(iNstar):
        for ipc in range(Nb):
            rPSF[ic]+=PCs[ipc]*icoeff[ic][ipc]
            #rPSF[ic]+=PCs[ipc]*coeffs[ic][ipc][0]
    return(rPSF)



def interp_cubic(image,target_Npix,dx,dy,osam):
    print(image.shape,target_Npix,dx,dy,osam)
    newf1=np.zeros(image.shape,dtype='float')
    Ng0=image.shape[0];Nge=target_Npix*osam
    dpix=((Ng0-Nge)/2);edx=dx*osam;edy=dy*osam
    outI=np.zeros((target_Npix,target_Npix),dtype='float')
    #print(dpix,edx,edy)
    for ic in range(Ng0):
        f1 = interpolate.interp1d(range(Ng0), image[ic] ,kind='cubic')
        newx=np.arange(0,Nge,1)+dpix-edx
        newf1[ic][0:Nge]=f1(newx)
    for ic in range(Nge):
        f2 = interpolate.interp1d(range(Ng0), newf1[:,ic], kind='cubic')
        newy=np.arange(0,Nge,1)+dpix-edy
        newf1[0:Nge,ic]=f2(newy)
    
    for i in range(target_Npix):
        for j in range(target_Npix):
            for it in range(osam):
                ra=int(i*osam+it);
                for jt in range(osam):
                    rb=int(j*osam+jt)
                    outI[i][j]+=newf1[ra][rb]
    return(outI)

            



def Hyb_interpolate(stamp_name,ofname,dirs,PCdirs,det,
    filte,PCs,pos,coeff,webcoeff,imodel,getNobj,center):
    print("in intepolate",stamp_name,ofname,dirs,PCdirs,det,filte)
    PCs=np.array(PCs);pos=np.array(pos);coeff=np.array(coeff)
    webcoeff=np.array(webcoeff);imodel=np.array(imodel);center=np.array(center)
    tmpname=stamp_name.split("/");l0=len(tmpname);l1=len(tmpname[l0-1]);
    osam=2
    dirs=dirs+"NIRCam/"+filte+"/"+det+"/"
    PCdirs=PCdirs+"NIRCam/"+filte+"/"+det+"/"
    tmpname=stamp_name.split("/");l0=len(tmpname);l1=len(tmpname[l0-1]);
    slicename=tmpname[l0-1][0:l1-len("_star.fits")];
    stamp_dirs=dirs+'/'+"catalogue/star_stamps/"
    cat_name=dirs+'/'+"catalogue/star/"+slicename+"_star.cat"
    mask_name=dirs+'/'+"catalogue/star_stamps/"+slicename+"_mask.fits"
    cat=np.loadtxt(cat_name);
    hdu=fits.open(ofname);OPD=hdu[0].header['DATE-BEG'];print("OPD:",OPD)
    PHOTUJA2=hdu['SCI'].header['PHOTUJA2'];XPOSURE=hdu['SCI'].header['XPOSURE']
    PHOTMJSR=hdu['SCI'].header['PHOTMJSR']
    PIXAR_A2=hdu['SCI'].header['PIXAR_A2']
    hdu=fits.open(mask_name);inmask=hdu[0].data;print("mask")
    hdu=fits.open(stamp_name);instar=hdu[0].data;print("stamp")
    Nobj=instar.shape[0]
    for ic in range(Nobj):
        instar[ic]=instar[ic].T*XPOSURE*PHOTMJSR
        inmask[ic]=inmask[ic].T
    Nobj=instar.shape[0];Ng0=instar.shape[1]
    spos=np.zeros((Nobj,2),dtype='float')
    wgt=np.zeros(instar.shape,dtype='float')
    for ic in range(Nobj):
        spos[ic][0]=cat[ic][1];spos[ic][1]=cat[ic][2]
        wgt[ic]=inmask[ic][:,:]
    print(PCdirs+slicename+"_mcPC.fits")
    write_mult_fits(PCdirs+slicename+"_mcPC.fits",PCs,pos,coeff,webcoeff);
    print(PCdirs+slicename+"_mcmodel.fits")
    write_mult_fits(PCdirs+slicename+"_mcmodel.fits",imodel,center);
    hdu=fits.open(PCdirs+'/'+slicename+"_mcPC.fits")
    hdu[0].header['filter']=filte
    hdu[0].header['detector']=detector[det].upper()
    hdu[0].header['data']=OPD
    hdu.writeto(PCdirs+'/'+slicename+"_mcPC.fits",overwrite=True)
    hdu.close()
    print(PCdirs+slicename+"_mcPC.fits")
    polyorder=2
    if getNobj>10:polyorder=3;print("1")
    rPSF=coeff2psf(spos,coeff,PCs,pos,polyorder);print("1",rPSF.shape)
    icoeff=polyfit2d(pos,webcoeff[:,0],webcoeff[:,1],pos,polyorder);print("2")
    eNg=rPSF.shape[1]
    Ng=int(eNg/osam)
    sNg=60;dpix=int((Ng0-sNg)/2);cut=int((imodel.shape[1]-eNg)/2);print("3")
    cent=np.zeros(2,dtype='float')
    mix_mod=np.zeros(rPSF.shape,dtype='float')
    chi1=np.zeros((getNobj,sNg,sNg))
    chi2=np.zeros((getNobj,sNg,sNg))
    residu1=np.zeros((getNobj,sNg,sNg))
    residu2=np.zeros((getNobj,sNg,sNg))
    fp=open(stamp_dirs+slicename+"_mcoshape.cat","w")
    fpw=open(stamp_dirs+slicename+"_mcwshape.cat","w")
    fpm=open(stamp_dirs+slicename+"_mcmshape.cat","w")
    for ic in range(getNobj):
        mean,sigma=gaus_estimate(instar[ic])
        print(ic,mean,sigma)
        instar[ic]-=mean
        instar[ic]*=wgt[ic]
        instar[ic]/=np.sum(instar[ic])
        #snrs=S2N(instar[ic]);
        tmp=wgtmap(instar[ic]);
        tmp_img=instar[ic][dpix:dpix+sNg,dpix:dpix+sNg]
        tmp_wgt=wgt[ic][dpix:dpix+sNg,dpix:dpix+sNg]#mask
        tmp_img*=tmp_wgt
        sums=np.sum(tmp_img);#tmp_img/=sums;#normalise
        print("sums",sums)
        error_map=tmp[dpix:dpix+sNg,dpix:dpix+sNg];error_map/=sums  #weight map
        tmp_wmap=1./(error_map*error_map)*tmp_wgt;
        cent[0]=center[ic][0]-dpix;cent[1]=center[ic][1]-dpix;
        print(cent[0],cent[1])
        dx=cent[0]-(sNg/2.+0.5)+1;dy=cent[1]-(sNg/2.+0.5)+1
        tmp_mod=interp_cubic(imodel[ic],sNg,dy,dx,osam)
        frac=(1.-np.sum(imodel[ic,cut:cut+eNg,cut:cut+eNg]))
        #rPSF=rPSF/np.sum(rPSF)*frac
        #tmp_mod*=tmp_wgt;tmp_mod/=np.sum(tmp_mod)
        for i in range(eNg):
            for j in range(eNg):
                mix_mod[ic][i][j]=0
                mix_mod[ic][i][j]=(rPSF[ic][i][j]+imodel[ic][i+cut][j+cut]*icoeff[ic])#mixture model
                if mix_mod[ic][i][j]<0 : mix_mod[ic][i][j]=0
        for i in range(eNg):
            for j in range(eNg):
                if mix_mod[ic][i][j]==0 :
                    mix_mod[ic][i][j]=(mix_mod[ic][i-1][j-1]+mix_mod[ic][i-1][j]+mix_mod[ic][i-1][j+1]+
                                       mix_mod[ic][i][j-1]+mix_mod[ic][i][j+1]+
                                       mix_mod[ic][i+1][j-1]+mix_mod[ic][i+1][j]+mix_mod[ic][i+1][j+1])/8.
        #tmp_mix_mod=model2pos(mix_mod[ic],sNg,cent,osam)
        tmp_mix_mod=interp_cubic(mix_mod[ic],sNg,dy,dx,osam)
        print(np.sum(tmp_img),np.sum(tmp_mix_mod))
        tmp_mix_mod*=tmp_wgt;tmp_mix_mod/=np.sum(tmp_mix_mod)
        residu1[ic]=(tmp_img-tmp_mod).T
        residu2[ic]=(tmp_img-tmp_mix_mod).T
        oe1,oe2,or2=size(tmp_img,cent,4.);
        we1,we2,wr2=size(tmp_mod,cent,4.);
        me1,me2,mr2=size(tmp_mix_mod,cent,4.);
        print(oe1,oe2,or2)
        print(me1,me2,mr2)
        fpw.writelines(str(we1)+'\t'+str(we2)+'\t'+str(wr2)+'\t'+str(spos[ic][0])+'\t'+str(spos[ic][1])+'\n')
        fp.writelines(str(oe1)+'\t'+str(oe2)+'\t'+str(or2)+'\t'+str(spos[ic][0])+'\t'+str(spos[ic][1])+'\n')
        fpm.writelines(str(me1)+'\t'+str(me2)+'\t'+str(mr2)+'\t'+str(spos[ic][0])+'\t'+str(spos[ic][1])+'\n')
        chi1[ic]=residu1[ic]*residu1[ic]*tmp_wmap
        chi2[ic]=residu2[ic]*residu2[ic]*tmp_wmap
    fp.close();fpw.close();fpm.close()
    write_fits(stamp_dirs+slicename+"_webchi.fits",chi1)
    write_fits(stamp_dirs+slicename+"_mixchi.fits",chi2)
    write_fits(stamp_dirs+slicename+"_webresidu.fits",residu1)
    write_fits(stamp_dirs+slicename+"_mixresidu.fits",residu2)  
    oshape=np.loadtxt(stamp_dirs+slicename+"_mcoshape.cat")
    wshape=np.loadtxt(stamp_dirs+slicename+"_mcmshape.cat")
    de1=oshape[:,0]-wshape[:,0];de2=oshape[:,1]-wshape[:,1];dr=(oshape[:,2]**0.5-wshape[:,2]**0.5)
    print(np.mean(oshape[:,0]),np.std(oshape[:,0]),np.mean(oshape[:,1]),np.std(oshape[:,1]),np.mean(oshape[:,2]),np.std(oshape[:,2]))
    print(np.mean(de1),np.std(de1),np.mean(de2),np.std(de2),np.median(dr),np.std(dr),np.median((oshape[:,2]-wshape[:,2])/oshape[:,2]))
    

