import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os
import ctypes
import sys
from scipy import interpolate
import matplotlib as mpl
import galsim
cmap = mpl.cm.hsv
import webbpsf
loaddir=os.getcwd()
libpsffit=ctypes.CDLL(loaddir+'/libpsffit.so')
from C_tools_lib import centriod_psf,Interp_bicubic,BInterp_bicubic

def write_fits(fitsname,data,head1=None,head2=None):
    wdata=data;
    if head1==None :
        #print("head1 none",fitsname)
        hdu = fits.PrimaryHDU(wdata)
        hdul = fits.HDUList([hdu])
        hdul.writeto(fitsname,overwrite=True)
    if head1!=None and head2==None:
        #print("head2 none",fitsname)
        hdu = fits.PrimaryHDU(wdata,header=head1)
        hdul = fits.HDUList([hdu])
        hdul.writeto(fitsname,overwrite=True)
    if head1!=None and head2!=None:
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

def get_filename(file_dir,append):
    fname=[];
    for list in os.listdir(file_dir):
        if os.path.isfile(os.path.join(file_dir,list)):
            if list[(len(list)-len(append)):len(list)]==append:
                fname.append(list);

    return fname



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

def psf_fit(stars,spos,npc,gain,method,wgt=None,osam=1,SNRs=100.):
    '''
    in star image stamps: stars[Nstar][Ng][Ng]
    in star image positions: spos[Nstar][0]:x;spos[Nstar][1]:y
    number of PCs
    gain of CCD: gain
    method of reconstructions
    '''
    Nstar=stars.shape[0];
    Ng0=stars.shape[1];
    dpix=19
    Nh=Ng0/2
    Ng=Ng0-dpix;
    print(Nstar,Ng0,npc,gain,method)
    fpoint=ctypes.POINTER(ctypes.c_float)
    ipoin1=ctypes.POINTER(ctypes.c_int)
    cint=ctypes.c_int
    cfloat=ctypes.c_float
    libpsffit.psf_fit.argtypes=\
    [fpoint,fpoint,cint,cint,cint,cint,fpoint,fpoint,cfloat,cint,ipoin1,cint,cfloat]

    star=((Nstar)*Ng0*Ng0*ctypes.c_float)()
    mask=((Nstar)*Ng0*Ng0*ctypes.c_float)()
    PCs=(npc*Ng*osam*Ng*osam*ctypes.c_float)()
    rePCs=np.zeros((npc,Ng*osam,Ng*osam),dtype=np.float)
    starpos=(Nstar*2*ctypes.c_float)()
    
    Coeff=(Nstar*npc*2*ctypes.c_float)()
    
    Nobj=(1*ctypes.c_int)()
    means=(Nstar*ctypes.c_float)()
    

    for ic in range(Nstar):
        k=ic*2
        starpos[k]=spos[ic][0]
        starpos[k+1]=spos[ic][1]
        for i in range(Ng0):
            for j in range(Ng0):
                k=ic*Ng0*Ng0+i*Ng0+j
                star[k]=stars[ic][i][j]

    if wgt is None:
      mask.all=1
    else:
      for ic in range(Nstar):
            for i in range(Ng0):
                for j in range(Ng0):
                    k=ic*Ng0*Ng0+i*Ng0+j
                    mask[k]=wgt[ic][i][j]

    libpsffit.psf_fit(star,starpos,npc,Nstar,Ng0,Ng,PCs,Coeff,gain,method,Nobj,osam,SNRs,mask)
    #the side need to be changed
    reNobj=Nobj[0]
    restarpos=np.zeros((reNobj,2),dtype=np.float)
    reCoeff=np.zeros((reNobj,npc,2),dtype=np.float)
    slecstar=np.zeros((reNobj,Ng0,Ng0),dtype=np.float)
    for ic in range(npc):
        for i in range(Ng*osam):
            for j in range(Ng*osam):
                k=ic*Ng*osam*Ng*osam+i*Ng*osam+j
                rePCs[ic][i][j]=PCs[k]
    
    for ic in range(reNobj):
        for i in range(npc):
            k=ic*npc*2+i*2
            reCoeff[ic][i][0]=Coeff[k]
            reCoeff[ic][i][1]=Coeff[k+1]
            #print(reCoeff[ic][i][0],)
        #print("\n")
        for i in range(2):
            k=ic*2
            restarpos[ic][0]=starpos[k];
            restarpos[ic][1]=starpos[k+1]
        for i in range(Ng0):
            for j in range(Ng0):
                slecstar[ic][i][j]=star[ic*Ng0*Ng0+i*Ng0+j]
    
    #print(reCoeff)
    return rePCs,restarpos,reCoeff,reNobj,slecstar

def web_psf_fit(stars,spos,npc,gain,method,wgt=None,imodel=None,osam=1,SNRs=100.):
    '''
    in star image stamps: stars[Nstar][Ng][Ng]
    in star image positions: spos[Nstar][0]:x;spos[Nstar][1]:y
    number of PCs
    gain of CCD: gain
    method of reconstructions
    '''
    print("osam=%d"%(osam))
    Nstar=stars.shape[0];
    Ng0=stars.shape[1];
    dpix=4
    Nh=Ng0/2
    Ng=Ng0-dpix;
    Ngm=int((imodel.shape[1])/osam);print("Ngm",Ngm)
    print(Nstar,Ng0,npc,gain,method)
    fpoint=ctypes.POINTER(ctypes.c_float)
    ipoin1=ctypes.POINTER(ctypes.c_int)
    cint=ctypes.c_int
    cfloat=ctypes.c_float
    libpsffit.web_psf_fit.argtypes=\
    [fpoint,fpoint,cint,cint,cint,cint,fpoint,fpoint,cfloat,cint,ipoin1,cint,cfloat,fpoint,fpoint,fpoint]

    star=((Nstar)*Ng0*Ng0*ctypes.c_float)()
    mask=((Nstar)*Ng0*Ng0*ctypes.c_float)()
    PCs=(npc*Ng*osam*Ng*osam*ctypes.c_float)()
    rePCs=np.zeros((npc,Ng*osam,Ng*osam),dtype='float')
    starpos=(Nstar*2*ctypes.c_float)()
    webmodel=((Nstar)*Ngm*osam*Ngm*osam*ctypes.c_float)()
    wcoeff=(Nstar*2*ctypes.c_float)()
    
    Coeff=(Nstar*npc*2*ctypes.c_float)()
    
    Nobj=(1*ctypes.c_int)()
    

    for ic in range(Nstar):
        k=ic*2
        starpos[k]=spos[ic][0]
        starpos[k+1]=spos[ic][1]
        for i in range(Ng0):
            for j in range(Ng0):
                k=ic*Ng0*Ng0+i*Ng0+j
                star[k]=stars[ic][i][j]


    if wgt is None:
        mask.all=1
    else:
        for ic in range(Nstar):
            for i in range(Ng0):
                for j in range(Ng0):
                    k=ic*Ng0*Ng0+i*Ng0+j
                    mask[k]=wgt[ic][i][j]

    if imodel is None:
        webmodel.all=0.
    else:
        for ic in range(Nstar):
            for i in range(Ngm*osam):
                for j in range(Ngm*osam):
                    k=ic*Ngm*osam*Ngm*osam+i*Ngm*osam+j
                    webmodel[k]=imodel[ic][i][j]

    print("start delivery")
    libpsffit.web_psf_fit(star,starpos,npc,Nstar,Ng0,Ng,PCs,Coeff,gain,method,Nobj,osam,SNRs,webmodel,mask,wcoeff)

    #the side need to be changed
    reNobj=Nobj[0]
    restarpos=np.zeros((reNobj,2),dtype=np.float)
    reCoeff=np.zeros((reNobj,npc,2),dtype=np.float)
    slecstar=np.zeros((reNobj,Ng0,Ng0),dtype=np.float)
    slecmask=np.zeros((reNobj,Ng0,Ng0),dtype=np.float)
    slecmodel=np.zeros((reNobj,Ngm*osam,Ngm*osam),dtype=np.float)
    webcoeff=np.zeros((reNobj,2),dtype=np.float)
    for ic in range(npc):
        for i in range(Ng*osam):
            for j in range(Ng*osam):
                k=ic*Ng*osam*Ng*osam+i*Ng*osam+j
                rePCs[ic][i][j]=PCs[k]
    
    for ic in range(reNobj):
        for i in range(npc):
            k=ic*npc*2+i*2
            reCoeff[ic][i][0]=Coeff[k]
            reCoeff[ic][i][1]=Coeff[k+1]
            #print(reCoeff[ic][i][0],)
        #print("\n")
        for i in range(2):
            k=ic*2
            restarpos[ic][0]=starpos[k];
            restarpos[ic][1]=starpos[k+1]
            webcoeff[ic][0]=wcoeff[k];
            webcoeff[ic][1]=wcoeff[k+1];
        for i in range(Ng0):
            for j in range(Ng0):
                slecstar[ic][i][j]=star[ic*Ng0*Ng0+i*Ng0+j]
                slecmask[ic][i][j]=mask[ic*Ng0*Ng0+i*Ng0+j]
        for i in range(Ngm*osam):
            for j in range(Ngm*osam):
                k=ic*Ngm*osam*Ngm*osam+i*Ngm*osam+j
                slecmodel[ic][i][j]=webmodel[k]
    
    #print(reCoeff)
    return rePCs,restarpos,reCoeff,reNobj,slecstar,webcoeff,slecmask,slecmodel



def psf_rec(psfpos,PCs,spos,Coeff,polyorder,method):
  '''
  interpolated PSF position: psfpos
  have 'poly', 'pixel' and 'krig'method
  '''
  Nobj=psfpos.shape[0];print("Nobj=",Nobj)
  Ng=PCs.shape[1]
  npc=PCs.shape[0]
  Nstar=Coeff.shape[0]
  fpoint=ctypes.POINTER(ctypes.c_float)
  ipoin1=ctypes.POINTER(ctypes.c_int)
  libpsffit.psf_rec.argtypes=[ctypes.c_int,ctypes.c_int,ctypes.c_int,fpoint,fpoint,\
  fpoint,ctypes.c_int,fpoint,ctypes.c_int,fpoint,ctypes.c_int]

  inpsfpos=(Nobj*2*ctypes.c_float)()
  inpsfimg=(Nobj*Ng*Ng*ctypes.c_float)()
  psfimg=np.zeros((Nobj,Ng,Ng),dtype=np.float)
  inPCs=(npc*Ng*Ng*ctypes.c_float)()
  inspos=(Nstar*2*ctypes.c_float)()
  inCoeff=(Nstar*npc*2*ctypes.c_float)()
  cmethod=ctypes.c_int
  if method=='poly':
    cmethod=1;print("poly")
  if method=='pixel':
    cmethod=2;print("pixel")
  if method=='krig':
    cmethod=3;print("krig")
  #print("shape:",spos.shape)
  for ic in range(Nstar):
    k=ic*2;
    inspos[k]=spos[ic][0]
    inspos[k+1]=spos[ic][1]
    #print(spos[ic][0],spos[ic][1],ic+1)
    for i in range(npc):
      k=ic*npc*2+i*2
      inCoeff[k]=Coeff[ic][i][0]
      inCoeff[k+1]=Coeff[ic][i][1]

  for ic in range(Nobj):
    k=ic*2
    inpsfpos[k]=psfpos[ic][0]
    inpsfpos[k+1]=psfpos[ic][1]
    for i in range(Ng):
      for j in range(Ng):
        k=ic*Ng*Ng+i*Ng+j
        inpsfimg[k]=psfimg[ic][i][j]
  for ic in range(npc):
    for i in range(Ng):
      for j in range(Ng):
        k=ic*Ng*Ng+i*Ng+j
        inPCs[k]=PCs[ic][i][j]

  print("start psf_rec")
  libpsffit.psf_rec(npc,polyorder,Nstar,inspos,inCoeff,inPCs,Ng,inpsfimg,Nobj,inpsfpos,cmethod)

  for ic in range(Nobj):
    for i in range(Ng):
      for j in range(Ng):
        k=ic*Ng*Ng+i*Ng+j
        psfimg[ic][i][j]=inpsfimg[k]
  
  return psfimg


def web_psf_rec(psfpos,PCs,spos,Coeff,polyorder,method):
  '''
  interpolated PSF position: psfpos
  have 'poly', 'pixel' and 'krig'method
  psfpos:required position of PSF
  '''
  rl=len(psfpos.shape);
  if rl==1 :
    Nobj=1;print("Nobj=",Nobj)
    Ng=PCs.shape[1]
    npc=PCs.shape[0]
    Nstar=Coeff.shape[0]
    fpoint=ctypes.POINTER(ctypes.c_float)
    ipoin1=ctypes.POINTER(ctypes.c_int)
    libpsffit.web_psf_rec.argtypes=[ctypes.c_int,ctypes.c_int,ctypes.c_int,fpoint,fpoint,\
    fpoint,ctypes.c_int,fpoint,ctypes.c_int,fpoint,ctypes.c_int]

    inpsfpos=(Nobj*2*ctypes.c_float)()
    inpsfimg=(Nobj*Ng*Ng*ctypes.c_float)()
    psfimg=np.zeros((Nobj,Ng,Ng),dtype=np.float)
    inPCs=(npc*Ng*Ng*ctypes.c_float)()
    inspos=(Nstar*2*ctypes.c_float)()
    inCoeff=(Nstar*npc*2*ctypes.c_float)()
    cmethod=ctypes.c_int
    if method=='poly':
        cmethod=1;print("poly")
    if method=='pixel':
        cmethod=2;print("pixel")
    if method=='krig':
        cmethod=3;print("krig")
    #print("shape:",spos.shape)
    for ic in range(Nstar):
        k=ic*2;
        inspos[k]=spos[ic][0]
        inspos[k+1]=spos[ic][1]
        #print(spos[ic][0],spos[ic][1],ic+1)
        for i in range(npc):
            k=ic*npc*2+i*2
            inCoeff[k]=Coeff[ic][i][0]
            inCoeff[k+1]=Coeff[ic][i][1]

    for ic in range(Nobj):
      k=ic*2
      inpsfpos[k]=psfpos[0]
      inpsfpos[k+1]=psfpos[1]
      for i in range(Ng):
        for j in range(Ng):
          k=ic*Ng*Ng+i*Ng+j
          inpsfimg[k]=psfimg[ic][i][j]
    for ic in range(npc):
        for i in range(Ng):
            for j in range(Ng):
                k=ic*Ng*Ng+i*Ng+j
                inPCs[k]=PCs[ic][i][j]

    print("start psf_rec")
    libpsffit.web_psf_rec(npc,polyorder,Nstar,inspos,inCoeff,inPCs,Ng,inpsfimg,Nobj,inpsfpos,cmethod)

    for ic in range(Nobj):
        for i in range(Ng):
            for j in range(Ng):
                k=ic*Ng*Ng+i*Ng+j
                psfimg[ic][i][j]=inpsfimg[k]
  
    return psfimg
  else :
    Nobj=psfpos.shape[0];print("Nobj=",Nobj)
    Ng=PCs.shape[1]
    npc=PCs.shape[0]
    Nstar=Coeff.shape[0]
    fpoint=ctypes.POINTER(ctypes.c_float)
    ipoin1=ctypes.POINTER(ctypes.c_int)
    libpsffit.web_psf_rec.argtypes=[ctypes.c_int,ctypes.c_int,ctypes.c_int,fpoint,fpoint,\
    fpoint,ctypes.c_int,fpoint,ctypes.c_int,fpoint,ctypes.c_int]

    inpsfpos=(Nobj*2*ctypes.c_float)()
    inpsfimg=(Nobj*Ng*Ng*ctypes.c_float)()
    psfimg=np.zeros((Nobj,Ng,Ng),dtype=np.float)
    inPCs=(npc*Ng*Ng*ctypes.c_float)()
    inspos=(Nstar*2*ctypes.c_float)()
    inCoeff=(Nstar*npc*2*ctypes.c_float)()
    cmethod=ctypes.c_int
    if method=='poly':
        cmethod=1;print("poly")
    if method=='pixel':
        cmethod=2;print("pixel")
    if method=='krig':
        cmethod=3;print("krig")
    #print("shape:",spos.shape)
    for ic in range(Nstar):
        k=ic*2;
        inspos[k]=spos[ic][0]
        inspos[k+1]=spos[ic][1]
        #print(spos[ic][0],spos[ic][1],ic+1)
        for i in range(npc):
            k=ic*npc*2+i*2
            inCoeff[k]=Coeff[ic][i][0]
            inCoeff[k+1]=Coeff[ic][i][1]

    for ic in range(Nobj):
        k=ic*2
        inpsfpos[k]=psfpos[ic][0]
        inpsfpos[k+1]=psfpos[ic][1]
        for i in range(Ng):
            for j in range(Ng):
                k=ic*Ng*Ng+i*Ng+j
                inpsfimg[k]=psfimg[ic][i][j]
    for ic in range(npc):
        for i in range(Ng):
            for j in range(Ng):
                k=ic*Ng*Ng+i*Ng+j
                inPCs[k]=PCs[ic][i][j]

    print("start psf_rec")
    libpsffit.web_psf_rec(npc,polyorder,Nstar,inspos,inCoeff,inPCs,Ng,inpsfimg,Nobj,inpsfpos,cmethod)

    for ic in range(Nobj):
        for i in range(Ng):
            for j in range(Ng):
                k=ic*Ng*Ng+i*Ng+j
                psfimg[ic][i][j]=inpsfimg[k]
  
    return psfimg


def scipy_interp(ispos,PCs,spos,coeff,method='bspline'):
    print("scipy_interp")
    from scipy import interpolate
    '''
    interpolated PSF positions by using scipy interpolate
    method:'bspline','rbf','2d','Nearest','Linear'
    '''
    npc=PCs.shape[0];Nstar=spos.shape[0];iNstar=ispos.shape[0];Ng=PCs.shape[1]
    icoeff=np.zeros((npc,iNstar),dtype='float')
    x=spos[:,0];y=spos[:,1];
    #fig=plt.figure(figsize=(40,30))
    if method=='bspline':
        print('bspline')
        for i in range(npc):
            z=coeff[:,i,0];W=1./coeff[:,i,1]
            tck=interpolate.bisplrep(x,y,z)
            for k in range(iNstar):
                icoeff[i][k]=interpolate.bisplev(ispos[k][0],ispos[k][1],tck)
    if method=='rbf':
        print('rbf')
        for i in range(npc):
            z=coeff[:,i,0];W=1./coeff[:,i,1]
            rbfi=interpolate.Rbf(x,y,z,function='thin_plate')
            for k in range(iNstar):
                icoeff[i][k]=rbfi(ispos[k][0],ispos[k][1])
    if method=='2d':
        print('2d')
        for i in range(npc):
            z=coeff[:,i,0];W=1./coeff[:,i,1]
            Spline=interpolate.interp2d(x,y,z,kind='linear')
            for k in range(iNstar):
                icoeff[i][k]=Spline(ispos[k][0],ispos[k][1])
    if method=='Nearest':
        print('Nearest')
        for i in range(npc):
            z=coeff[:,i,0];W=1./coeff[:,i,1]
            interp =interpolate.NearestNDInterpolator(spos, z)
            for k in range(iNstar):
                icoeff[i][k]=interp(ispos[k][0],ispos[k][1])
    if method=='Linear':
        print('Linear')
        for i in range(npc):
            z=coeff[:,i,0];W=1./coeff[:,i,1]
            Lint=interpolate.LinearNDInterpolator(list(zip(x, y)), z)
            for k in range(iNstar):
                icoeff[i][k]=Lint(ispos[k][0],ispos[k][1])
    '''for i in range(npc):
        z=coeff[:,i,0];W=1./coeff[:,i,1]
        #S=interpolate.SmoothBivariateSpline(x,y,z)
        #tck=interpolate.bisplrep(x,y,z)
        #Spline=interpolate.interp2d(x,y,z,kind='linear')
        rbfi=interpolate.Rbf(x,y,z,function='inverse')
        #Lint=interpolate.LinearNDInterpolator(spos, z)
        for k in range(iNstar):
            #icoeff[i][k]=S(ispos[k][0],ispos[k][1])
            #icoeff[i][k]=interpolate.bisplev(ispos[k][0],ispos[k][1],tck)
            #icoeff[i][k]=Spline(ispos[k][0],ispos[k][1])
            icoeff[i][k]=rbfi(ispos[k][0],ispos[k][1])
            #icoeff[i][k]=Lint(ispos[k][0],ispos[k][1])
        vmin=np.min(z);vmin1=np.min(icoeff[i,:]);mins=min(vmin,vmin1)
        vmax=np.max(z);vmax1=np.max(icoeff[i,:]);maxs=max(vmax,vmax1)
        zp=(z-mins)/(maxs-mins);zp1=(icoeff[i,:]-mins)/(maxs-mins)
        plt.subplot(3,4,i+1)
        for k in range(Nstar):
            plt.plot(x[k],y[k],'.',c=cmap(zp[k]))
        for k in range(iNstar):
            plt.plot(ispos[k][0]+100,ispos[k][1]+100,'*',c=cmap(zp1[k]))
    #plt.savefig('surface.pdf')
    '''
    restar=np.zeros((iNstar,Ng,Ng))
    for ic in range(iNstar):
        for i in range(npc):
            restar[ic]=restar[ic]+PCs[i]*icoeff[i][ic]
    return(restar)


def Shepard(ispos,PCs,spos,coeff,indx):
    print("Shepard")
    Nobj=ispos.shape[0];Nstar=spos.shape[0]
    npc=PCs.shape[0];Ng=PCs.shape[1]
    icoeff=np.zeros((npc,Nobj),dtype='float')
    W=np.zeros((Nstar),dtype='float')
    print(W.shape[0],coeff.shape[0],icoeff.shape[0])
    for ic in range(Nobj):
        for k in range(Nstar):
            denominator=np.sqrt((ispos[ic][0]-spos[k][0])**2+(ispos[ic][1]-spos[k][1])**2)
            if denominator==0:W[k]=10.**10
            else:
                W[k]=(np.sqrt((ispos[ic][0]-spos[k][0])**2+(ispos[ic][1]-spos[k][1])**2))**(-indx)
        for l in range(npc):
            icoeff[l][ic]=np.dot(W,coeff[:,l,0])/np.sum(W)

    restar=np.zeros((Nobj,Ng,Ng))
    for ic in range(Nobj):
        for i in range(npc):
            restar[ic]=restar[ic]+PCs[i]*icoeff[i][ic]
    return(restar)

def psfex_rec(ispos,pex,Ng):
    print("psfex")
    import galsim
    import galsim.des as gdes
    iNstar=ispos.shape[0]
    des_psfex = gdes.DES_PSFEx(pex)
    rpsf=np.zeros((iNstar,Ng,Ng),dtype='float')
    for ic in range(iNstar):
        image_pos=galsim.PositionD(ispos[ic][0], ispos[ic][1])
        psf=des_psfex.getPSF(image_pos)
        rpsf[ic]=psf.image.array[:,:]
    return(rpsf)



def estimate_para(stars,center):
    Nobj=stars.shape[0]
    Ng=stars.shape[1]
    dpoint=ctypes.POINTER(ctypes.c_double)
    instar=(Nobj*Ng*Ng*ctypes.c_double)()
    incent=(Nobj*2*ctypes.c_double)()
    inpara=(5*ctypes.c_double)()
    libpsffit.estimate_parapy.argtypes=[dpoint,dpoint,ctypes.c_int,ctypes.c_int,\
    ctypes.c_int,dpoint]
    for ic in range(Nobj):
        #cx,cy=centriod_psf(stars[ic])
        #print(cx,cy)
        for i in range(Ng):
            for j in range(Ng):
                instar[ic*Ng*Ng+i*Ng+j]=stars[ic][i][j]
        incent[ic*2]=center[ic][0];incent[ic*2+1]=center[ic][1]
    print("start estimate_para")
    libpsffit.estimate_parapy(instar,incent,Nobj,Ng,5,inpara)
    sigma=inpara[0]
    rd=inpara[1]
    beta=inpara[2]
    return(sigma,rd,beta)
    
def star_shape(stars,pos,fname,rg=3.):
    Nobj=stars.shape[0];
    Ng=stars.shape[1]
    center=np.zeros((Nobj,2),dtype='float')
    rg=0;
    for ic in range(Nobj):
        center[ic][0],center[ic][1],sig=centriod_psf(stars[ic])
        rg=rg+sig
    rg,rd,beta=estimate_para(stars,center)
    #rg=rg/Nobj*3.5;
    print("rg=%.5f"%(rg));rg=3.5
    fp=open(fname,"w")
    for ic in range(Nobj):
        e1,e2,R2=size(stars[ic],center[ic],rg)
        fp.writelines(str(R2)+'\t'+str(e1)+'\t'+str(e2)+'\t'+str(pos[ic][0])+'\t'+str(pos[ic][1])+'\n')
    fp.close()

def starcat2psfcat(starcat):
    ''' we assume that the psf = star/np.sum(star)
    parameters:
    ---------------------
    starcat [3D array], e.g. N star, starcat.shape = (N, 42, 42) 
    returns:
    ---------------------
    psfcat [3D array], starcat.shape = psfcat.shape
    '''
    star_shape = starcat.shape
    x, y, z = star_shape
    star_sum = np.sum(np.sum(starcat, axis=1), axis=1)
    starcat_reshape = np.ravel(starcat).reshape(x, y*z).T
    psfcat = (np.ravel((starcat_reshape/star_sum).T)).reshape(star_shape)
    return psfcat

def fft_conv(image,psf):
    '''
    convlove image with a psf,usually image size should be larger than psf
    '''
    size=image.shape
    core=np.zeros(size,dtype='float')
    nx=image.shape[0];ny=image.shape[1]
    px=psf.shape[0];py=psf.shape[1]
    ra=int((nx-px)/2);rb=int((ny-py)/2)
    core[ra-1:ra+px-1,rb-1:rb+py-1]=psf
    Fg=np.fft.fft2(image)
    Fp=np.fft.fft2(core)
    Fg*=Fp
    for i in range(Fp.real.shape[0]):
        for j in range(Fp.real.shape[1]):
            l=(i+j+2)%2
            if (l==1) :Fg.real[i][j]*=-1;Fg.imag[i][j]*=-1;
    tmp=np.fft.ifft2(Fg)
    cimage=tmp.real
    return(cimage)

def int2cent(image,Target_Ng):
  Ng0=image.shape[0]
  cent=np.zeros(2,dtype='float')
  cent[0],cent[1],sigma=centriod_psf(image)
  Target=Interp_bicubic(image,cent,Target_Ng,Target_Ng)
  return(Target)


def gaus_estimate(image):
    """
    gaus_estimate(image)
    return(mean,sigma)
    """
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




def model2pos(image,Target_Ng,pos,osam):
    cent=np.zeros(2,dtype='float')
    cent[0]=pos[0]*osam+osam/2.
    cent[1]=pos[1]*osam+osam/2.
    align=BInterp_bicubic(image,cent,Target_Ng*osam,Target_Ng*osam)
    Target=np.zeros((Target_Ng,Target_Ng),dtype='float')
    for i in range(Target_Ng):
        for j in range(Target_Ng):
            for l in range(osam):
                ra=i*osam+l
                for k in range(osam):
                    rb=j*osam+k
                    Target[i][j]=Target[i][j]+align[ra][rb]
    return(Target)

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
    #print(image.shape,target_Npix,dx,dy,osam)
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



def get_webPSF(PCfile,
              psfpos,
              Target_Ng,
              osam=1,
              sign="M"):
    '''
    the keywords of function get_webPSF are:
    PCfile,       the constructed PC file
    psfpos,         the required psf pos which you want, pos can be a 2 dimentional array,such as psfpos[xpos][ypos]
    Target_Ng,      the required pixel size of psf stamp
    osam=1,         the oversampleing factor of the PSF, osam can be 1 or 2
    '''
    print("in get_webPSF")
    if sign=="M":rNstar=6
    else : rNstar=100
    import webbpsf
    tmpname=PCfile.split("/");l0=len(tmpname);l1=len(tmpname[l0-1]);
    slicename=tmpname[l0-1][0:l1-len(".fits")]
    try:
        hdu=fits.open(PCfile);
    except:
        print("file: "+PCfile+" open failed, please check the file name!")
        sys.exit()
    current_dir=os.getcwd()
    print("here")
    Nobj=psfpos.shape[0];psfcat=np.zeros((Nobj,2));print("Nobj:",Nobj)
    for ic in range(Nobj):
        psfcat[ic][0]=psfpos[ic][0];psfcat[ic][1]=psfpos[ic][1];
    rPSF=np.zeros((Nobj,Target_Ng*osam,Target_Ng*osam),dtype='float')
    Tmodel=np.zeros((Nobj,Target_Ng*osam,Target_Ng*osam),dtype='float')
    if os.path.exists(PCfile)==False :
        Nstar=0
        OPDs=hdu[0].header['DATE']
        Filter=hdu[0].header['FILTER']
        Detect=hdu[0].header['DETECTOR']
        Osam=hdu[0].header['OSAM']
        webp=webbpsf.NIRCam();
        webp.filter=Filter;
        webp.detector=Detect
        webp.load_wss_opd_by_date(OPDs)
        for ic in range(Nobj):
            print(ic)
            webp.detector_position=(psfcat[ic][0],psfcat[ic][1])
            psf = webp.calc_psf(fov_pixels=Target_Ng,oversample=osam)
            psf=psf[2].data.T;psf/=np.sum(psf)
            Tmodel[ic]=psf
    else:
        # use mixture psf model
        OPDs=hdu[0].header['DATE']
        Filter=hdu[0].header['FILTER']
        Detect=hdu[0].header['DETECTOR']
        Osam=hdu[0].header['OSAM']
        PCs=hdu[0].data
        npc=PCs.shape[0]
        spos=hdu[1].data
        coeff=hdu[2].data
        webcoeff=hdu[3].data
        Nstar=spos.shape[0]
        print(Nstar,rNstar)
        if Nstar<rNstar:
            print("too few star, use webbpsf")
            webp=webbpsf.NIRCam();
            webp.filter=Filter;
            webp.detector=Detect
            webp.load_wss_opd_by_date(OPDs)
            for ic in range(Nobj):
                print(ic)
                webp.detector_position=(psfcat[ic][0],psfcat[ic][1])
                psf = webp.calc_psf(fov_pixels=Target_Ng,oversample=osam)
                psf=psf[2].data.T;psf/=np.sum(psf)
                Tmodel[ic]=psf
        else:
            print("hybird psf model")
            coeff=hdu[2].data;print(psfpos.shape)
            if Nstar>10:polyorder=3
            if Nstar>15:polyorder=4
            if Nstar>21:polyorder=5
            rPSF=coeff2psf(spos,coeff,PCs,psfpos,degrees=polyorder);print("rPSF",rPSF.shape)
            eNg=rPSF.shape[1];print("eNg=",eNg)
            sNg=eNg;dNg=int((eNg-sNg)/2)
            cut_rPSF=np.zeros(rPSF.shape,dtype='float')
            for ic in range(Nobj):
                for i in range(eNg):
                    ra=i-(eNg/2-0.5)
                    for j in range(eNg):
                        rb=j-(eNg/2-0.5)
                        r2=ra**2+rb**2
                        r=r2**0.5
                        if r2<(sNg/2)**2:
                            cut_rPSF[ic][i][j]=1
                            cut_rPSF[ic][i][j]=1./(1+np.exp((3.*(r-9))))
                cut_rPSF[ic]*=rPSF[ic]
            icoeff=polyfit2d(spos,webcoeff[:,0],webcoeff[:,1],psfpos,polyorder=3)
            #webcoeff[:,0]
            print(icoeff)
            if Target_Ng<eNg/Osam:
                print("Target_Ng:%d\t%d"%(Target_Ng,eNg/2));
                cent=np.zeros(2,dtype='float');
                webp=webbpsf.NIRCam();
                webp.filter=Filter;
                webp.detector=Detect
                webp.load_wss_opd_by_date(OPDs)
                dpix=int((eNg-Target_Ng*Osam)/2);print("dpix:",dpix)
                for ic in range(Nobj):
                    webp.detector_position=(psfcat[ic][0],psfcat[ic][1])
                    psf = webp.calc_psf(fov_pixels=eNg/Osam,oversample=Osam)
                    psf=psf[2].data.T;psf/=np.sum(psf)
                    #psf*=0
                    '''plt.imshow((np.fabs(cut_rPSF[ic,:,:]))**0.1);plt.show()
                    plt.imshow((np.fabs(psf[:,:]*icoeff[ic]))**0.1);plt.show()'''
                    tmppsf=psf[:,:]*icoeff[ic]+cut_rPSF[ic,:,:]
                    #plt.imshow((tmppsf)**0.1);plt.show()
                    for i in range(eNg):
                        for j in range(eNg):
                            if  tmppsf[i][j]<0 : #padding the nagetive value
                                tmppsf[i][j]=0
                    for i in range(eNg):
                        for j in range(eNg):
                            if  tmppsf[i][j]==0 and i-1>0 and i+1<eNg and j-1>0 and j+1<eNg: #padding the nagetive value
                                tmppsf[i][j]=(tmppsf[i-1][j-1]+tmppsf[i-1][j]+tmppsf[i-1][j+1]+
                                    tmppsf[i][j-1]+tmppsf[i][j+1]+
                                    tmppsf[i+1][j-1]+tmppsf[i+1][j]+tmppsf[i+1][j+1])/8.
                    tmppsf/=np.sum(tmppsf)
                    print("osam=%d"%(osam))
                    #write_fits("fits/test.fits",tmppsf)
                    if osam==2 and Osam==2:
                        tmp=tmppsf[dpix:dpix+Target_Ng*osam,dpix:dpix+Target_Ng*osam]
                        Tmodel[ic]=tmp/np.sum(tmp)
                    if osam==1 and Osam==2:#merger pixels
                        print("case 1")
                        Tmodel[ic]*=0
                        tmp=tmppsf[dpix:dpix+Target_Ng*2,dpix:dpix+Target_Ng*2]
                        for i in range(Target_Ng):
                            for j in range(Target_Ng):
                                for l in range(2):
                                    ra=i*2+l;
                                    for m in range(2):
                                        rb=j*2+m
                                        Tmodel[ic][i][j]+=tmp[ra][rb]
                    if osam==1 and Osam==1:#merger pixels
                        print("case 1")
                        Tmodel[ic]*=0
                        tmp=tmppsf[dpix:dpix+Target_Ng,dpix:dpix+Target_Ng]
                        Tmodel[ic]=tmp
                        Tmodel[ic]/=np.sum(Tmodel[ic])
                    if osam==2 and Osam==1:
                        print("the oversample factor of PC is 1, re generate PC with osam=2")
                        sys.exit()
                    #write_fits("fits/test.fits",Tmodel)
            else:  #Target_Ng>=eNg/osam
                dpix=int((Target_Ng*Osam-eNg)/2);print("case 2")
                webp=webbpsf.NIRCam();
                webp.filter=Filter;
                webp.detector=Detect
                webp.load_wss_opd_by_date(OPDs)
                for ic in range(Nobj):
                    webp.detector_position=(psfcat[ic][0],psfcat[ic][1]) 
                    psf = webp.calc_psf(fov_pixels=Target_Ng,oversample=Osam)
                    psf=psf[2].data.T;psf/=np.sum(psf)
                    psf/=np.sum(psf[dpix:dpix+eNg,dpix:dpix+eNg]) #calibrate the flux ratio
                    psf*=icoeff[ic]
                    #psf*=0
                    psf[dpix:dpix+eNg,dpix:dpix+eNg]+=cut_rPSF[ic,:,:]
                    #psf*=0;psf[dpix:dpix+eNg,dpix:dpix+eNg]=cut_rPSF[ic,:,:]
                    tmp=psf
                    #plt.imshow(((tmp))**0.1);plt.show()
                    for i in range(eNg):
                        for j in range(eNg):
                            if  tmp[i+dpix][j+dpix]<0 : #padding the nagetive value
                                tmp[i+dpix][j+dpix]=0
                    #plt.imshow(((tmp))**0.1);plt.show()
                    for i in range(Target_Ng*Osam):
                        for j in range(Target_Ng*Osam):
                            if  tmp[i][j]==0 and i-1>0 and i+1<Target_Ng*Osam and j-1>0 and j+1<Target_Ng*Osam: #padding the nagetive value
                                tmp[i][j]=(tmp[i-1][j-1]+tmp[i-1][j]+tmp[i-1][j+1]+tmp[i][j-1]+tmp[i][j+1]+tmp[i+1][j-1]+tmp[i+1][j]+tmp[i+1][j+1])/8.
                    tmp/=np.sum(tmp)
                    #plt.imshow(((tmp))**0.1);plt.show()
                    if osam==2 and Osam==2:
                        Tmodel[ic,:,:]=tmp[:,:]
                    if osam==1 and Osam==2:#merger pixels
                        Tmodel[ic]*=0
                        for i in range(Target_Ng):
                            for j in range(Target_Ng):
                                for l in range(2):
                                    ra=i*2+l;
                                    for m in range(2):
                                        rb=j*2+m
                                        Tmodel[ic][i][j]+=tmp[ra][rb]
                    if osam==1 and Osam==1:
                        Tmodel[ic]*=0
                        Tmodel[ic]=tmp
                    if osam==2 and Osam==1:
                        print("the oversample factor of PC is 1, re generate PC with osam=2")
                        sys.exit()
                    Tmodel[ic]/=np.sum(Tmodel[ic]) 
    return(Tmodel)


















