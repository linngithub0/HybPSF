import numpy as np
import ctypes
import os
from astropy.io import fits
#pulic declaretions
pwd=os.getcwd()
libpsffit=ctypes.CDLL(pwd+'/libpsffit.so')
fpoint=ctypes.POINTER(ctypes.c_float)
dpoint=ctypes.POINTER(ctypes.c_double)
ipoint=ctypes.POINTER(ctypes.c_int)



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



def centriod_psf(psf):
	'''
	intput the psf arry with Ng,Ng,
	and return the estimated pixel center of the profile :cx,cy
	based on fast shape estimate algorithm 
	'''
	Ng1=psf.shape[0]
	Ng2=psf.shape[1];
	in_image=(Ng1*Ng2*ctypes.c_float)(*list(np.ravel(psf)))
	cent=(3*ctypes.c_double)()
	libpsffit.centriod_psflib.argtypes=[fpoint,ctypes.c_int,ctypes.c_int,dpoint]
	libpsffit.centriod_psflib(in_image,Ng1,Ng2,cent)
	return(cent[0],cent[1],cent[2])
def centriod_galaxy(galaxy):
	'''
	intput the galaxy arry with Ng,Ng,
	and return the estimated pixel center of the profile :cx,cy
	based on fast shape estimate algorithm 
	'''
	Ng1=galaxy.shape[0]
	Ng2=galaxy.shape[1];
	in_image=(Ng1*Ng2*ctypes.c_float)(*list(np.ravel(galaxy)))
	cent=(2*ctypes.c_double)()
	libpsffit.centriod_galaxylib.argtypes=[fpoint,ctypes.c_int,ctypes.c_int,dpoint]
	libpsffit.centriod_galaxylib(in_image,Ng1,Ng2,cent)
	return(cent[0],cent[1])

def centriod_burst(galaxy):
	'''
	intput the galaxy arry with Ng,Ng,
	and return the estimated pixel center of the profile :cx,cy
	based on fast shape estimate algorithm 
	'''
	Ng1=galaxy.shape[0]
	Ng2=galaxy.shape[1];
	in_image=(Ng1*Ng2*ctypes.c_float)(*list(np.ravel(galaxy)))
	cent=(2*ctypes.c_double)()
	libpsffit.centriod_burstlib.argtypes=[fpoint,ctypes.c_int,ctypes.c_int,dpoint]
	libpsffit.centriod_burstlib(in_image,Ng1,Ng2,cent)
	return(cent[0],cent[1])


def Interp_bicubic(image,center,cutx,cuty):
	'''
	using the pixel center of image and cubic spline to interpolate the image
	to the center of the pixel frame
	'''
	Ng1=image.shape[0]
	Ng2=image.shape[1]
	interpolated=(cutx*cutx*ctypes.c_float)()
	cent=(2*ctypes.c_double)(*list(np.ravel(center)))
	in_image=(Ng1*Ng2*ctypes.c_float)(*list(np.ravel(image)))
	libpsffit.Interp_bicubiclib.argtypes=[fpoint,ctypes.c_int,ctypes.c_int,\
	ctypes.c_int,ctypes.c_int,dpoint,fpoint]
	libpsffit.Interp_bicubiclib(in_image,Ng1,Ng2,cutx,cuty,cent,interpolated)
	out=np.array(list(interpolated)).reshape(cutx,cuty)[0:cutx]
	return(out)

def BInterp_bicubic(image,center,cutx,cuty):
	'''
	using the pixel center of image and cubic spline to interpolate the image
	to the center of the pixel frame
	'''
	Ng1=image.shape[0]
	Ng2=image.shape[1]
	interpolated=(cutx*cutx*ctypes.c_float)()
	cent=(2*ctypes.c_double)(*list(np.ravel(center)))
	in_image=(Ng1*Ng2*ctypes.c_float)(*list(np.ravel(image)))
	libpsffit.BInterp_bicubiclib.argtypes=[fpoint,ctypes.c_int,ctypes.c_int,\
	ctypes.c_int,ctypes.c_int,dpoint,fpoint]
	libpsffit.BInterp_bicubiclib(in_image,Ng1,Ng2,cutx,cuty,cent,interpolated)
	out=np.array(list(interpolated)).reshape(cutx,cuty)[0:cutx]
	return(out)


def Noise_padding(image1,image2,mean):
    Ngx=image1.shape[0];Ngy=image1.shape[1]
    val=np.empty((0,1),float);
    r2min=(min(Ngx,Ngy)/2.)**2;
    for i in range(Ngx):
        for j in range(Ngy):
            x=i-Ngx/2.+0.5;
            y=j-Ngy/2.+0.5
            c=image2[i][j]
            r2=x*x+y*y
            if r2>r2min and c!=0:
                val=np.append(val,[[c],],axis=0)
    lens=val.shape[0];#print(lens)
    for i in range(Ngx):
        for j in range(Ngy):
            indx=np.random.randint(low=0,high=lens-1)
            if image1[i][j]==0: image1[i][j]=val[indx]+mean
    return(image1)


def Deblend(image):
    Ngx=image.shape[0];Ngy=image.shape[1]
    mean,sigma=gaus_estimate(image);print(mean,sigma)
    image-=mean
    in_sigma=(1*ctypes.c_float)()
    blendings=(3*ctypes.c_int)()
    in_image=(Ngx*Ngy*ctypes.c_float)(*list(np.ravel(image)))
    out_image=(Ngx*Ngy*ctypes.c_float)(*list(np.ravel(image)))
    in_sigma=sigma
    libpsffit.deblending.argtypes=[fpoint,fpoint,ctypes.c_int,ctypes.c_int,ctypes.c_float,ipoint]
    libpsffit.deblending(in_image,out_image,Ngx,Ngy,in_sigma,blendings)
    out=np.array(list(out_image)).reshape(Ngx,Ngy)[0:Ngx]
    ngroup=blendings[0];nsub1=blendings[1];nsub2=blendings[2]
    #mean,sigma=gaus_estimate(out);
    out1=Noise_padding(out,image,mean)
    return(out1,ngroup,nsub1,nsub2)


























