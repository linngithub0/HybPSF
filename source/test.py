import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os
import time
import webbpsf

from psf_fit import write_fits,fft_conv,get_filename
from assist import cat_header,gains
from config import file_name,cata_path,Camra,Filter,Detector,PC_path,\
npc,diagnose,osam

from psf_fit import star_shape,size,write_fits,get_filename,\
web_psf_fit,web_psf_rec,int2cent,model2pos,wgtmap,gaus_estimate,write_mult_fits
from C_tools_lib import centriod_psf



tmpname=file_name.split("/");l0=len(tmpname);l1=len(tmpname[l0-1])
slicename=tmpname[l0-1][0:l1-5];
fname=tmpname[l0-1];
if l0 == 1:
    file_path=libdir=os.getcwd();
else:
    file_path=file_name[0:(len(file_name)-len(fname)-1)]

catdirs=wdirs=cata_path+'/'
fdirs=file_path
#generate catalogue by using sextractor

hdu=fits.open(fdirs+fname) 
try:
    Instr=hdu[0].header['INSTRUME'];print(Instr)
    Detec=hdu[0].header['DETECTOR'];print(Detec)
    Filt=hdu[0].header['FILTER'];print(Filt)
except:
    from config import Instr,Detec,Filt
data=hdu[1].data                     #get fits data
wgt=np.ones(data.shape,dtype='float')
coldcat=wdirs+"cold.cat"
segname=wdirs+slicename+"_seg.fits"
coldcmd="sex "+file_name+" -CATALOG_NAME "+coldcat+" -DETECT_MINAREA 120 -DETECT_THRESH 3.2 \
-DEBLEND_NTHRESH 32 -DEBLEND_MINCONT 0.04 -BACK_SIZE 400 -BACK_FILTERSIZE 5 -BACKPHOTO_TYPE LOCAL \
-CHECKIMAGE_TYPE SEGMENTATION -CHECKIMAGE_NAME "+segname
print(coldcmd)
os.system(coldcmd)
hotcat=wdirs+"hot.cat"
hotcmd="sex "+file_name+" -CATALOG_NAME "+hotcat+" -DETECT_MINAREA 18 -DETECT_THRESH 1.2 \
-DEBLEND_NTHRESH 32 -DEBLEND_MINCONT 0.065 -BACK_SIZE 100 -BACK_FILTERSIZE 3 -BACKPHOTO_TYPE LOCAL"
print(hotcmd)
os.system(hotcmd)

try:  
    #start erosion
    erosion_core=np.ones((40,40))
    hdu=fits.open(segname)
    segmap=hdu[0].data
    erosion=fft_conv(segmap,erosion_core)
    cold_cat=np.loadtxt(coldcat);print("cold:",cold_cat.shape)
    hot_cat=np.loadtxt(hotcat);print("hot:",hot_cat.shape)
    Nobj=hot_cat.shape[0]
    hcopy=np.zeros(hot_cat.shape,dtype='float')
    k=0
    for it in range(Nobj):
        x=int((hot_cat[it][1]));y=int((hot_cat[it][2]))
        if erosion[x][y]<0 :
            hcopy[k]=hot_cat[it]
            k=k+1
    cold_cat=np.append(cold_cat,hcopy[0:k,:],axis=0)
    indx=np.where(cold_cat[:,5]>0)
    cold_cat=cold_cat[indx[0]]#select only have flux larger than 0
    print("total slected source number:",cold_cat.shape[0])
    fp=open(wdirs+slicename+".cat","w")
    for line in cat_header:
        fp.writelines(line)
    for it in range(cold_cat.shape[0]):
        fp.writelines(str(it+1)+'\t'+str(cold_cat[it][1])+'\t'\
            +str(cold_cat[it][2])+'\t'+str(cold_cat[it][3])+'\t'\
            +str(cold_cat[it][4])+'\t'+str(cold_cat[it][5])+'\t'\
            +str(cold_cat[it][6])+'\t'+str(cold_cat[it][7])+'\t'\
            +str(cold_cat[it][8])+'\n')
    fp.close()
    #r50=cold_cat[:,4];fwhm=cold_cat[:,3];mag=-2.5*np.log10(cold_cat[:,5])
    #fig=plt.figure()
    #plt.plot(mag,r50,'.',c='blue',alpha=0.3)
    #plt.plot(mag,fwhm,'.',c='red',alpha=0.5)
    #plt.ylim(-2,30)
    #plt.savefig(wdirs+slicename+'_cat.pdf')
except:
    print("file open failed, please check this file %s"%(fname))
    pass


staridf=("./staridf/cut "+catdirs+" "+slicename+".cat "
    +fdirs+" "+tmpname[l0-1]+" "+catdirs+" -8. -3. "+str(gains[Detec]))
print(staridf)
print("try find star in : %s"%(catdirs+slicename+".cat"))
os.system(staridf)
cat=catdirs+slicename+".cat"
cold_cat=np.loadtxt(cat)
r50=cold_cat[:,4];fwhm=cold_cat[:,3];mag=-2.5*np.log10(cold_cat[:,5])
rcat=catdirs+"catalogue/star/"+slicename+"_star.cat"
print(rcat)
cat_fits=catdirs+"catalogue/star_stamps/"+slicename+"_star.fits"
hdu=fits.open(cat_fits);print(hdu.info())
hdr=hdu[0].header
hdr['INSTRUME']=Instr
hdr['FILTER']=Filt
hdr['DETECTOR']=Detec
hdu.writeto(cat_fits,overwrite=True)
hdu.close()
try:
    rcat=np.loadtxt(rcat);
    rfwhm=rcat[:,3];rmag=-2.5*np.log10(rcat[:,5]);rr50=rcat[:,4];
    fig=plt.figure()
    plt.plot(mag,r50,'.',c='blue',alpha=0.3)
    plt.plot(mag,fwhm,'.',c='red',alpha=0.3)
    plt.plot(rmag,rfwhm,'^',c='red',alpha=0.3)
    plt.plot(rmag,rr50,'*',c='blue',alpha=0.3)
    plt.plot(rmag[0],rr50[0],'*',c='blue',alpha=0.3,label='R50')
    plt.plot(rmag[0],rfwhm[0],'*',c='blue',alpha=0.3,label='FWHM')
    plt.ylim(-2,30)
    plt.legend()
    plt.title(slicename+".fits"+"(%d stars)"%(rcat.shape[0]))
    plt.savefig(catdirs+slicename+'_cat.pdf')
    print("the star stamps are now stored in path:"+wdirs+"catalogue/star_stamps/"+
        ", and with name: "+slicename+"_star.fits")
    print("the star catalogue is now generated at path:"+wdirs+"catalogue/star/"+
        ", and with name: "+slicename+"_star.cat")
except:
    print(catdirs+slicename+".fits does not have enough stars")
    pass



stamp_dirs=wdirs+"catalogue/star_stamps/"
stamp_name=stamp_dirs+slicename+"_star.fits"
mask_name=stamp_dirs+slicename+"_mask.fits"
wgt_name=wdirs+slicename+"_wgt.fits"
cat_name=wdirs+"catalogue/star/"+slicename+"_star.cat"

try: 
    cat=np.loadtxt(cat_name);
    hdu=fits.open(stamp_name);OPD=hdu[0].header['DATE'];print("OPD:",OPD)
    instar=hdu[0].data;print("instar shape:",instar.shape)
    Nobj=instar.shape[0]  
    hdu=fits.open(mask_name)
    inmask=hdu[0].data;
except :
    print("file:"+stamp_name+" does not exist!, PSF model are constructed by webbpsf")
    Nobj=0
    pass

if Nobj<npc+1:
    print("Nobj=%d, no enough star find in this catalogue"%(Nobj));
    pass
else:
    #use mixture model......
    Ng0=instar.shape[1]
    Nh=instar.shape[1]/2
    for ic in range(Nobj):
        instar[ic]=instar[ic].T
        inmask[ic]=inmask[ic].T
    spos=np.zeros((Nobj,2),dtype='float')
    for ic in range(cat.shape[0]):
        spos[ic][0]=cat[ic][1];spos[ic][1]=cat[ic][2]
        rx=cat[ic][1]+0.5;ry=cat[ic][2]+0.5
        xl=int(rx-Nh);xh=int(rx+Nh);yl=int(ry-Nh);yh=int(ry+Nh)
    imodel=np.zeros((Nobj,(Ng0-1)*osam,(Ng0-1)*osam),dtype='float')
    if Instr=='NIRCAM':
            nrc = webbpsf.NIRCam()
            nrc.filter =Filt
            nrc.detector=Detec
            nrc.load_wss_opd_by_date(OPD)
    else:
        print("only NIRCAM can be used at present code!")
        exit()
    for ic in range(Nobj):
        print(slicename+" mixpsf imodel",ic+1)
        nrc.detector_position=(spos[ic][0],spos[ic][1])
        psf = nrc.calc_psf(fov_pixels=Ng0-1,oversample=2)
        imodel[ic]=(psf[0].data/np.sum(psf[0].data)).T
    wgt=inmask
    PCs,pos,coeff,getNobj,slecstar=web_psf_fit(
        instar,
        spos,
        npc,
        gain=float(gains[Detec]),
        method=1,
        wgt=wgt,
        imodel=imodel,
        osam=2,
        SNRs=0.0000001)
    write_mult_fits(wdirs+slicename+"_PC.fits",PCs,pos,coeff)
    hdu=fits.open(wdirs+slicename+"_PC.fits")
    hdr=hdu[0].header
    hdr['INSTRUME']=Instr
    hdr['FILTER']=Filt
    hdr['DETECTOR']=Detec
    hdu.writeto(cat_fits,overwrite=True)
    hdu.close()
    print("The PC data are now stored at path: "+wdirs+"with name: "+slicename+"_PC.fits")
    if diagnose==True:
        rPSF=web_psf_rec(spos,PCs,pos,coeff,polyorder=1,method='poly')
        eNg=rPSF.shape[1]
        Ng=int(eNg/osam)
        for ic in range(getNobj):
            for i in range(eNg):
                for j in range(eNg):
                    if np.fabs(rPSF[ic][i][j])<10**(-4):rPSF[ic][i][j]=0 #drop nagetive value
        sNg=int(Ng0/2);dpix=int(Ng0/4);cut=int(((Ng0-1)*osam-eNg)/2)
        cent=np.zeros(2,dtype='float')
        mix_mod=np.zeros(rPSF.shape,dtype='float')
        chi1=np.zeros((getNobj,sNg,sNg))
        chi2=np.zeros((getNobj,sNg,sNg))
        residu1=np.zeros((getNobj,sNg,sNg))
        residu2=np.zeros((getNobj,sNg,sNg))
        wmodel=np.zeros((getNobj,sNg,sNg))
        mmodel=np.zeros((getNobj,sNg,sNg))
        fp=open(wdirs+slicename+"_oshape.cat","w")
        fpw=open(wdirs+slicename+"_wshape.cat","w")
        fpm=open(wdirs+slicename+"_mshape.cat","w")
        for ic in range(getNobj):
            mean,sigma=gaus_estimate(instar[ic])
            instar[ic]-=mean
            tmp=wgtmap(instar[ic],gain=1)
            tmp_img=instar[ic][dpix:dpix+sNg,dpix:dpix+sNg]
            tmp_wgt=wgt[ic][dpix:dpix+sNg,dpix:dpix+sNg]#mask
            tmp_img*=tmp_wgt
            sums=np.sum(tmp_img);tmp_img/=sums #normalise
            error_map=tmp[dpix:dpix+sNg,dpix:dpix+sNg];error_map/=sums  #weight map
            tmp_wmap=1./(error_map*error_map)*tmp_wgt;
            cent[0],cent[1],sigma=centriod_psf(tmp_img);print("cx=%f,cy=%f"%(cent[0],cent[1]))
            tmp_mod=model2pos(imodel[ic],sNg,cent,osam)
            tmp_mod*=tmp_wgt;tmp_mod/=np.sum(tmp_mod)
            for i in range(eNg):
                for j in range(eNg):
                    mix_mod[ic][i][j]=0
                    mix_mod[ic][i][j]+=rPSF[ic][i][j]+imodel[ic][i+cut][j+cut]#mixture model
            for i in range(eNg):
                for j in range(eNg):
                    if mix_mod[ic][i][j]<0 :
                         mix_mod[ic][i][j]=(mix_mod[ic][i-1][j-1]+mix_mod[ic][i-1][j]+mix_mod[ic][i-1][j+1]+
                            mix_mod[ic][i][j-1]+mix_mod[ic][i][j+1]+
                            mix_mod[ic][i+1][j-1]+mix_mod[ic][i+1][j]+mix_mod[ic][i+1][j+1])/8.
            tmp_mix_mod=model2pos(mix_mod[ic],sNg,cent,osam)
            tmp_mix_mod*=tmp_wgt;tmp_mix_mod/=np.sum(tmp_mix_mod)
            residu1[ic]=tmp_img-tmp_mod;wmodel[ic]=tmp_mod;
            residu2[ic]=tmp_img-tmp_mix_mod;mmodel[ic]=tmp_mix_mod
            oe1,oe2,or2=size(tmp_img,cent,4.);fp.writelines(str(oe1)+'\t'+str(oe2)+'\t'+str(or2)+'\n')
            we1,we2,wr2=size(tmp_mod,cent,4.);fpw.writelines(str(we1)+'\t'+str(we2)+'\t'+str(wr2)+'\n')
            me1,me2,mr2=size(tmp_mix_mod,cent,4.);fpm.writelines(str(me1)+'\t'+str(me2)+'\t'+str(mr2)+'\n')
            chi1[ic]=residu1[ic]*residu1[ic]*tmp_wmap
            chi2[ic]=residu2[ic]*residu2[ic]*tmp_wmap
        fp.close();fpw.close();fpm.close()
        write_fits(wdirs+slicename+"_webchi.fits",chi1)
        write_fits(wdirs+slicename+"_mixchi.fits",chi2)
        write_fits(wdirs+slicename+"_webresidu.fits",residu1)
        write_fits(wdirs+slicename+"_mixresidu.fits",residu2)
        write_fits(wdirs+slicename+"_wmodel.fits",wmodel)
        write_fits(wdirs+slicename+"_mmodel.fits",mmodel)
        print("the diagnose file are write to "+wdirs)






