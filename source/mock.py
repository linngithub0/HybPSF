import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os
import ctypes
import sys
import matplotlib as mpl
cmap = mpl.cm.hsv
import sys
import webbpsf


from psf_fit import psf_fit,psf_rec,scipy_interp,Shepard,star_shape,\
size,write_fits,get_filename,psfex_rec

nrc = webbpsf.NIRCam()  #get the instrument
print("load NIRCam down")
nrc.filter='F115W'      #declare the filter
nrc.detector='NRCA3'    #declare the detector
Nobj=100
nx=ny=2040
spos=np.zeros((Nobj,2),dtype='float')
instar=np.zeros((Nobj,Ng0,Ng0),dtype='float')
imodel=np.zeros((Nobj,Ng0*osam,Ng0*osam),dtype='float')
spos[:,0]=np.random.uniform(0,nx,Nobj)
spos[:,1]=np.random.uniform(0,ny,Nobj)


print("start creat examples")
f1=open("fits/pos.cat","w")
for ic in range(Nobj):
	print(ic)
	nrc.detector_position=(spos[ic][0],spos[ic][1])
	psf=nrc.calc_psf(fov_pixels=Ng0,oversample=2)
	instar[ic]=psf[1].data
	imodel[ic]=psf[0].data
	f1.writelines(str(spos[ic][0])+"\t"+str(spos[ic][1])+"\n")
f1.close()

write_fits("fits/instar.fits",instar)
write_fits("fits/imodel.fits",imodel)