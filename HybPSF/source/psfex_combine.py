import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os
import ctypes
import sys
import matplotlib as mpl
cmap = mpl.cm.hsv

from psf_fit import write_fits,get_filename

#psfex
#convert multiple .psfex into one
num=int(sys.argv[1])
if num<10 :
	chip='0'+str(num)
else :
	chip=str(num)
dirs="/Volumes/SHIN/Nie/psfTest/TEST_150/MSC_00000"
wdirs="/Volumes/SHIN/Nie/psfTest/TEST_150/"
SNRs="300"
append=chip+"_raw.fits.cat"
count=0
chip_num=15
fnames=[];psfcata=[]
print(fnames)
for i in range(chip_num):
	print(i)
	if i<9 :
		fdir=dirs+"0"+str(i+1)+"/"
		fname=get_filename(fdir,append);
		for j in range(len(fname)):
			if fname[j][0]=='r' and fname[j][1]=='e':
				fnames.append(fdir+fname[j])
	else:
		fdir=dirs+str(i+1)+"/"
		fname=get_filename(fdir,append)
		for j in range(len(fname)):
			if fname[j][0]=='r' and fname[j][1]=='e':
				fnames.append(fdir+fname[j])

print(fnames)
hdul=fits.open(fnames[0]);print("len",len(fnames))
whdul=hdul.copy();print("count",hdul[2].header['NAXIS2'])
for i in range(len(fnames)-1):
	hdul=fits.open(fnames[i+1]);print("count",hdul[2].header['NAXIS2'])
	whdul[2].data=np.append(whdul[2].data,hdul[2].data)
	count=whdul[2].header['NAXIS2']

whdul[2].header=hdul[2].header
whdul[2].header['NAXIS2']=count
print(count)
wname="/Volumes/SHIN/Nie/psfTest/TEST_150/"+chip+".cat"
whdul.writeto(wname,overwrite=True)
Ng=27
import os
shell_cmd="psfex "+wname+" -SAMPLE_MINSN "+SNRs+" -PSF_DIR "+wdirs+SNRs+"/ -PSF_SIZE "+str(Ng)+","+str(Ng)
os.system(shell_cmd)