import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import ctypes
import webbpsf


from psf_fit import star_shape,size,write_fits,get_filename,\
web_psf_fit,web_psf_rec,int2cent,model2pos,wgtmap,gaus_estimate
from C_tools_lib import centriod_psf



osam=2
npc=3
gain=1
Cam="NIRCam/"
filte="F090W/"
obs="jw02736001001_02101_"
detect=["nrca1/","nrca2/","nrca3/","nrca4/","nrcb1/","nrcb2/","nrcb3/","nrcb4/"]
dirs="/home/lnie/data/JWST/"+Cam+filte
fdirs=dirs+detect[0]
stamp_dirs=dirs+detect[0]+"catalogue/star_stamps/"
get_name=get_filename(fdirs,"_cal.fits")
stamp_name=stamp_dirs+get_name[0]+".fits"
tmpname=get_name[0].split("/");l0=len(tmpname);l1=len(tmpname[l0-1]);slicename=tmpname[l0-1][0:l1-5]
wgt_name=fdirs+slicename+"_wgt.fits"
cat_name=dirs+detect[0]+"catalogue/star/"+get_name[0]+".cat"
print(stamp_name,wgt_name,cat_name)

cat=np.loadtxt(cat_name)
hdu=fits.open(stamp_name)
instar=hdu[0].data
Nobj=instar.shape[0]
Ng0=instar.shape[1]
Nh=instar.shape[1]/2
for ic in range(Nobj):
    instar[ic]=instar[ic].T
wgt=np.zeros(instar.shape,dtype='float')
spos=np.zeros((Nobj,2),dtype='float')

hdu=fits.open(wgt_name)
dwgt=hdu[0].data;dwgt=dwgt.T

k=0
for ic in range(cat.shape[0]):
    spos[ic][0]=cat[ic][1];spos[ic][1]=cat[ic][2]
    rx=cat[ic][1]+0.5;ry=cat[ic][2]+0.5
    xl=int(rx-Nh);xh=int(rx+Nh);yl=int(ry-Nh);yh=int(ry+Nh)
    wgt[k]=dwgt[xl:xh,yl:yh]
    k=k+1
imodel=np.zeros((Nobj,(Ng0-1)*osam,(Ng0-1)*osam),dtype='float')
nrc = webbpsf.NIRCam()
nrc.filter =filte[0:len(filte)-1]
tmp=detect[0].upper();
nrc.detector=tmp[0:len(filte)-1]
print(filte[0:len(filte)-1],tmp[0:len(filte)-1])
for ic in range(Nobj):
    print("imodel",ic+1)
    nrc.detector_position=(spos[ic][0],spos[ic][1])
    psf = nrc.calc_psf(fov_pixels=Ng0-1,oversample=2)
    imodel[ic]=(psf[0].data/np.sum(psf[0].data)).T
instars=instar
inspos=spos
inwgt=wgt
PCs,pos,coeff,getNobj,slecstar=web_psf_fit(
    instars,
    inspos,
    npc,
    gain,
    method=1,
    wgt=inwgt,
    imodel=imodel,
    osam=2,
    SNRs=0.0000001)


write_fits("fits/PCs.fits",PCs)


'''
#interpolate at the reconstruction position
rPSF=web_psf_rec(spos,PCs,spos,coeff,polyorder=1,method='poly')
#rPSF represent the PSFs
eNg=rPSF.shape[1]
Ng=int(eNg/osam)
sNg=int(Ng0/2);dpix=int(Ng0/4);cut=int(((Ng0-1)*osam-eNg)/2)
cent=np.zeros(2,dtype='float')
mix_mod=np.zeros(rPSF.shape,dtype='float')
tmp_img=np.zeros((getNobj,sNg,sNg))
tmp_wgt=np.zeros((getNobj,sNg,sNg))
tmp_wmap=np.zeros((getNobj,sNg,sNg))
error_map=np.zeros((getNobj,sNg,sNg))
tmp_mix_mod=np.zeros((getNobj,sNg,sNg))
tmp_mod=np.zeros((getNobj,sNg,sNg))
chi1=np.zeros((getNobj,sNg,sNg))
chi2=np.zeros((getNobj,sNg,sNg))
residu1=np.zeros((getNobj,sNg,sNg))
residu2=np.zeros((getNobj,sNg,sNg))
for ic in range(getNobj):
	for i in range(rPSF.shape[1]):
		for j in range(rPSF.shape[1]):
			if np.fabs(rPSF[ic][i][j])<10**(-4):
				rPSF[ic][i][j]=0
for ic in range(getNobj):
	mean,sigma=gaus_estimate(instars[ic])
	print("sigma=%f,mean=%f"%(sigma,mean))
	instars[ic]-=mean
	tmp=wgtmap(instars[ic],gain=1)
	tmp_img[ic]=instars[ic][dpix:dpix+sNg,dpix:dpix+sNg]#image
	tmp_wgt[ic]=inwgt[ic][dpix:dpix+sNg,dpix:dpix+sNg]#mask
	tmp_img[ic]*=tmp_wgt[ic]
	sums=np.sum(tmp_img[ic])
	tmp_img[ic]/=sums
	error_map[ic]=tmp[dpix:dpix+sNg,dpix:dpix+sNg]#weight map
	error_map[ic]/=sums
	#error_map[ic]*=tmp_wgt[ic]
	tmp_wmap[ic]=1./(error_map[ic]*error_map[ic])*tmp_wgt[ic]
	#tmp_wmap[ic]*=sums*sums
	cent[0],cent[1],sigma=centriod_psf(tmp_img[ic])
	tmp_mod[ic]=model2pos(imodel[ic],sNg,cent,osam)#webbpsf model
	tmp_mod[ic]*=tmp_wgt[ic]
	tmp_mod[ic]/=np.sum(tmp_mod[ic])
	for i in range(eNg):
		for j in range(eNg):
			mix_mod[ic][i][j]+=rPSF[ic][i][j]+imodel[ic][i+cut][j+cut]#mixture model
	tmp_mix_mod[ic]=model2pos(mix_mod[ic],sNg,cent,osam)
	tmp_mix_mod[ic]*=tmp_wgt[ic]
	tmp_mix_mod[ic]/=np.sum(tmp_mix_mod[ic])
	residu1[ic]=tmp_img[ic]-tmp_mod[ic]
	residu2[ic]=tmp_img[ic]-tmp_mix_mod[ic]
	print("residu1=%e,residu2=%e"%(np.sum(residu1[ic]),np.sum(residu2[ic])))
	chi1[ic]=residu1[ic]*residu1[ic]*tmp_wmap[ic]
	chi2[ic]=residu2[ic]*residu2[ic]*tmp_wmap[ic]
	print("chi1=%e,chi2=%e"%(np.sum(chi1[ic])/(sNg*sNg),np.sum(chi2[ic])/(sNg*sNg)))
	#tmp_img[ic]*=sums

write_fits("fits/1.fits",tmp_img)
write_fits("fits/2.fits",tmp_wgt)
write_fits("fits/3.fits",tmp_wmap)
write_fits("fits/4.fits",tmp_mod)
write_fits("fits/5.fits",tmp_mix_mod)
write_fits("fits/residu1.fits",residu1)
write_fits("fits/residu2.fits",residu2)
write_fits("fits/chi1.fits",chi1)
write_fits("fits/chi2.fits",chi2)
write_fits("fits/error.fits",error_map)




#interpolated test


rpos=spos[getNobj:Nobj,:]
rstar=instar[getNobj:Nobj,:,:]
#write_fits("fits/rstar.fits",rstar)
rwgt=wgt[getNobj:Nobj,:,:]
rPSF=web_psf_rec(rpos,PCs,spos,coeff,polyorder=1,method='poly')
Nobj=rPSF.shape[0]

mix_mod=np.zeros(rPSF.shape,dtype='float')
tmp_img=np.zeros((Nobj,sNg,sNg))
tmp_wgt=np.zeros((Nobj,sNg,sNg))
tmp_wmap=np.zeros((Nobj,sNg,sNg))
error_map=np.zeros((Nobj,sNg,sNg))
tmp_mix_mod=np.zeros((Nobj,sNg,sNg))
tmp_mod=np.zeros((Nobj,sNg,sNg))
chi1=np.zeros((Nobj,sNg,sNg))
chi2=np.zeros((Nobj,sNg,sNg))
residu1=np.zeros((Nobj,sNg,sNg))
residu2=np.zeros((Nobj,sNg,sNg))

for ic in range(Nobj):
	for i in range(rPSF.shape[1]):
		for j in range(rPSF.shape[1]):
			if np.fabs(rPSF[ic][i][j])<10**(-5):
				rPSF[ic][i][j]=0
for ic in range(Nobj):
	mean,sigma=gaus_estimate(rstar[ic])
	rstar[ic]-=mean
	tmp=wgtmap(rstar[ic],gain=1)
	tmp_img[ic]=rstar[ic][dpix:dpix+sNg,dpix:dpix+sNg]#image
	tmp_wgt[ic]=rwgt[ic][dpix:dpix+sNg,dpix:dpix+sNg]#mask
	tmp_img[ic]*=tmp_wgt[ic]
	sums=np.sum(tmp_img[ic])
	tmp_img[ic]/=sums
	error_map[ic]=tmp[dpix:dpix+sNg,dpix:dpix+sNg]#weight map
	error_map[ic]/=sums
	tmp_wmap[ic]=1./(error_map[ic]*error_map[ic])*tmp_wgt[ic]
	cent[0],cent[1],sigma=centriod_psf(tmp_img[ic])
	#construct models
	nrc.detector_position=(rpos[ic][0],rpos[ic][1])
	psf = nrc.calc_psf(fov_pixels=Ng0-1,oversample=2)
	nmodel=(psf[0].data/np.sum(psf[0].data)).T
	tmp_mod[ic]=model2pos(nmodel,sNg,cent,osam)#webbpsf model
	tmp_mod[ic]*=tmp_wgt[ic]
	tmp_mod[ic]/=np.sum(tmp_mod[ic])
	for i in range(eNg):
		for j in range(eNg):
			mix_mod[ic][i][j]+=rPSF[ic][i][j]+nmodel[i+cut][j+cut]#mixture model
	tmp_mix_mod[ic]=model2pos(mix_mod[ic],sNg,cent,osam)
	tmp_mix_mod[ic]*=tmp_wgt[ic]
	tmp_mix_mod[ic]/=np.sum(tmp_mix_mod[ic])
	residu1[ic]=tmp_img[ic]-tmp_mod[ic]
	residu2[ic]=tmp_img[ic]-tmp_mix_mod[ic]
	print("residu1=%e,residu2=%e"%(np.sum(residu1[ic]),np.sum(residu2[ic])))
	chi1[ic]=residu1[ic]*residu1[ic]*tmp_wmap[ic]
	chi2[ic]=residu2[ic]*residu2[ic]*tmp_wmap[ic]
	print("chi1=%e,chi2=%e"%(np.sum(chi1[ic])/(sNg*sNg),np.sum(chi2[ic])/(sNg*sNg)))



write_fits("fits/1.fits",tmp_img)
write_fits("fits/2.fits",tmp_wgt)
write_fits("fits/3.fits",tmp_wmap)
write_fits("fits/4.fits",tmp_mod)
write_fits("fits/5.fits",tmp_mix_mod)
write_fits("fits/iresidu1.fits",residu1)
write_fits("fits/iresidu2.fits",residu2)
write_fits("fits/ichi1.fits",chi1)
write_fits("fits/ichi2.fits",chi2)
write_fits("fits/ierror.fits",error_map)
'''











