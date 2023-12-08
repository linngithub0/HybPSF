from astropy.io import fits
from psf_fit import write_fits,fft_conv
import os
import numpy as np
import matplotlib.pyplot as plt
cat_header=["#   1 NUMBER                 Running object number\n",                                  
"#   2 X_IMAGE                Object position along x                                      [pixel]\n",
"#   3 Y_IMAGE                Object position along y                                      [pixel]\n",
"#   4 FWHM_IMAGE             FWHM assuming a gaussian core                                [pixel]\n",
"#   5 FLUX_RADIUS            Fraction-of-light radii                                      [pixel]\n",
"#   6 FLUX_AUTO              Flux within a Kron-like elliptical aperture                  [count]\n",
"#   7 FLUXERR_AUTO           RMS error for AUTO flux                                      [count]\n",
"#   8 FLAGS                  Extraction flags"                                          
"#   9 CLASS_STAR             S/G classifier output\n"]

fdir="/Users/linn/Documents/data/JWST/jw02736001001_02101_00004_nrca1/"

fname=fdir+"jw02736001001_02101_00004_nrca1_cal.fits"
tmp=fname.split('/');l0=len(tmp);l1=len(tmp[l0-1])
rename=tmp[l0-1]
hdu=fits.open(fname)
data=hdu[1].data
head1=hdu[0].header
head2=hdu[1].header
wname=fdir+rename


write_fits(wname,data,head1=head1,head2=head2)#extract only the sci data into wname
wgt=np.ones(data.shape,dtype='float')
for i in range(data.shape[0]):
	for j in range(data.shape[1]):
		if data[i][j]==0 or data[i][j]<0:
			wgt[i][j]=0
write_fits(fdir+rename[0:l1-5]+"_wgt.fits",wgt)

cmd="sex "+wname+" -CATALOG_NAME "+fdir+"cold.cat"\
+" -DETECT_MINAREA 120 -DETECT_THRESH 3.2 -DEBLEND_NTHRESH 32 -DEBLEND_MINCONT 0.04\
 -BACK_SIZE 400 -BACK_FILTERSIZE 5 -BACKPHOTO_TYPE LOCAL -CHECKIMAGE_TYPE SEGMENTATION -CHECKIMAGE_NAME "+\
 fdir+"segmap.fits"
os.system(cmd)
print(cmd)



cmd="sex "+wname+" -CATALOG_NAME "+fdir+"hot.cat"\
+" -DETECT_MINAREA 18 -DETECT_THRESH 1.2 -DEBLEND_NTHRESH 32 -DEBLEND_MINCONT 0.065\
 -BACK_SIZE 100 -BACK_FILTERSIZE 3 -BACKPHOTO_TYPE LOCAL"
os.system(cmd)
print(cmd)


erosion_core=np.ones((40,40))
hdu=fits.open(fdir+"segmap.fits")
segmap=hdu[0].data
print("start fft_conv")
erosion=fft_conv(segmap,erosion_core)
print("fft_conv down")
erosion=erosion.T
#write_fits(fdir+"erosion.fits",erosion)
cold_cat=np.loadtxt(fdir+"cold.cat");print(cold_cat.shape)
hot_cat=np.loadtxt(fdir+"hot.cat");print(hot_cat.shape)
Nobj=hot_cat.shape[0]
hcopy=np.zeros(hot_cat.shape,dtype='float')
k=0
for ic in range(Nobj):
	x=int((hot_cat[ic][1]));y=int((hot_cat[ic][2]))
	if erosion[x][y]<0 :
		hcopy[k]=hot_cat[ic]
		k=k+1
cold_cat=np.append(cold_cat,hcopy[0:k,:],axis=0)
print(k,cold_cat.shape[0])
indx=np.where(cold_cat[:,5]>0)
cold_cat=cold_cat[indx[0]]#select only have flux larger than 0
print(k,cold_cat.shape[0])

#cold_cat=hot_cat
fp=open(fdir+rename[:l1-5]+".cat","w")
for line in cat_header:
	fp.writelines(line)
for ic in range(cold_cat.shape[0]):
#for ic in range(10):
	fp.writelines(str(ic+1)+'\t'+str(cold_cat[ic][1])+'\t'\
		+str(cold_cat[ic][2])+'\t'+str(cold_cat[ic][3])+'\t'\
		+str(cold_cat[ic][4])+'\t'+str(cold_cat[ic][5])+'\t'\
		+str(cold_cat[ic][6])+'\t'+str(cold_cat[ic][7])+'\t'\
		+str(cold_cat[ic][8])+'\n')

fp.close()
r50=cold_cat[:,4];fwhm=cold_cat[:,3];mag=-2.5*np.log10(cold_cat[:,5])
plt.plot(mag,r50,'.',c='blue',alpha=0.3)
plt.plot(mag,fwhm,'.',c='red',alpha=0.5)
plt.ylim(-2,30)
plt.savefig(fdir+'cat.pdf')


#./cut <raw cata file path> <input cata name> <fits path> <input fitsname> <output path>
# [<left> <right>] [<weight fits>]

cmd="./staridf/cut "+fdir+" "+rename[:l1-5]+".cat "+fdir+" "+rename+" "+fdir+" -8. 2."
print(cmd)
os.system(cmd)



rcat=fdir+"catalogue/star/"+rename+".cat"
rcat=np.loadtxt(rcat)
rfwhm=rcat[:,3];rmag=-2.5*np.log10(rcat[:,5]);rr50=rcat[:,4];
plt.plot(mag,r50,'.',c='blue',alpha=0.3)
plt.plot(mag,fwhm,'.',c='red',alpha=0.3)
plt.plot(rmag,rfwhm,'^',c='red',alpha=0.3)
plt.plot(rmag,rr50,'*',c='blue',alpha=0.3)
plt.ylim(-2,30)
plt.savefig(fdir+'rcat.pdf')




