import numpy as np
import ctypes
#pulic declaretions
libpsffit=ctypes.CDLL('./libpsffit.so')
fpoint=ctypes.POINTER(ctypes.c_float)
dpoint=ctypes.POINTER(ctypes.c_double)


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
