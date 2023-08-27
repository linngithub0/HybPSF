import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os
import time
import sys
import webbpsf
from scipy.stats import sigmaclip
from astropy.stats import SigmaClip


from psf_fit import write_fits,fft_conv,get_filename,coeff2psf,interp_cubic,\
polyfit2d,S2N
from assist import cat_header,gains,detector
from config import file_name,cata_path,npc,diagnose,osam

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
fdirs=file_path+'/'
#generate catalogue by using sextractor

hdu=fits.open(fdirs+fname) 
OPD=hdu[0].header['DATE-BEG']
print("OPD:",OPD)
PHOTUJA2=hdu['SCI'].header['PHOTUJA2'];
XPOSURE=hdu['SCI'].header['XPOSURE']
PHOTMJSR=hdu['SCI'].header['PHOTMJSR']
PIXAR_A2=hdu['SCI'].header['PIXAR_A2']


try:
    Instr=hdu[0].header['INSTRUME'];print(Instr)
    Detec=hdu[0].header['DETECTOR'];print(Detec)
    '''
    if Detec=="NRCALONG" : Detec="NRCA5"
    if Detec=="NRCBLONG" : Detec="NRCB5"
    '''
    Filt=hdu[0].header['FILTER'];print(Filt)
except:
    from config import Instr,Detec,Filt
try:
    data=hdu['SCI'].data                     #get fits data
    header_append="[SCI]"
except:
    data=hdu[1].data
    header_append="\0"
wgt=np.ones(data.shape,dtype='float')
coldcat=wdirs+"cold.cat"
segname=wdirs+slicename+"_seg.fits"
coldcmd="source-extractor "+file_name+" -CATALOG_NAME "+coldcat+" -DETECT_MINAREA 120 -DETECT_THRESH 3.2 \
-DEBLEND_NTHRESH 32 -DEBLEND_MINCONT 0.04 -BACK_SIZE 400 -BACK_FILTERSIZE 5 -BACKPHOTO_TYPE LOCAL \
-CHECKIMAGE_TYPE SEGMENTATION -CHECKIMAGE_NAME "+segname
print(coldcmd)
os.system(coldcmd)
hotcat=wdirs+"hot.cat"
hotcmd="source-extractor "+file_name+" -CATALOG_NAME "+hotcat+" -DETECT_MINAREA 18 -DETECT_THRESH 1.2 \
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
    fwhm=cold_cat[:,3];mag=-2.5*np.log10(cold_cat[:,5]);r50=cold_cat[:,4];
    plt.plot(mag,fwhm,'.',c='red',alpha=0.3)
    plt.plot(mag,r50,'.',c='blue',alpha=0.3)
    plt.ylim(-2,30)
    plt.savefig(catdirs+slicename+'_oat.pdf')
except:
    print("file open failed, please check this file %s"%(fname))
    pass


staridf=("./cut "+catdirs+" "+slicename+".cat "+fdirs+" "+tmpname[l0-1]
         +" "+catdirs+" -7. -2. "+str(gains[Detec]))
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
hdu[0].header['INSTRUME']=Instr
hdu[0].header['FILTER']=Filt
hdu[0].header['DETECTOR']=Detec
hdu[0].header['DATE-BEG']=OPD
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
Nim=5

try: 
    cat=np.loadtxt(cat_name)
    hdu=fits.open(stamp_name)
    readstar=hdu[0].data*XPOSURE*PHOTMJSR;
    print("instar shape:",readstar.shape)
    Nobj=readstar.shape[0]  
    hdu=fits.open(mask_name)
    readmask=hdu[0].data
except :
    print("file:"+stamp_name+" does not exist!, PSF model are constructed by webbpsf")
    Nobj=0
    pass

cent=np.zeros(2,dtype='float')

if Nobj<npc+1:
    print("Nobj=%d, no enough star find in this catalogue"%(Nobj));
    pass
else:
    #use mixture model......
    Ng0=readstar.shape[1]
    Nh=readstar.shape[1]/2
    for ic in range(Nobj):
        readstar[ic]=readstar[ic].T
        readmask[ic]=readmask[ic].T
    readspos=np.zeros((Nobj,2),dtype='float')
    readwgt=np.zeros(readstar.shape,dtype='float')
    shapes=np.zeros((Nobj,2),dtype='float')
    shiftxy=[]
    for ic in range(Nobj):
        readspos[ic][0]=cat[ic][1];readspos[ic][1]=cat[ic][2]
        rx=cat[ic][1]+0.5;ry=cat[ic][2]+0.5
        xl=int(rx-Nh);xh=int(rx+Nh);yl=int(ry-Nh);yh=int(ry+Nh)
        readwgt[ic]=readmask[ic][:,:]
        mean,sigma=gaus_estimate(readstar[ic]);
        cent[0],cent[1],sigma1=centriod_psf(readstar[ic]-mean)
        oe1,oe2,or2=size(readstar[ic],cent,4.);
        mage=(oe1**2+oe2**2)**0.5
        shiftxy.append([np.fabs(cent[0]-Ng0/2+0.5),np.fabs(cent[1]-Ng0/2+0.5)])
        shapes[ic][0]=mage;shapes[ic][1]=or2
    cleane,low,high=sigmaclip(shapes[:,0],3.,3.)
    cleanr,low,high=sigmaclip(shapes[:,1],3.,3.)
    shiftxy=np.array(shiftxy)
    #clean the outliers
    print("sigmaclip clean Nobj=",cleane.shape[0])
    k=0
    for ic in range(Nobj):
        if (k <cleane.shape[0] and shapes[ic][0]==cleane[k]):
            mean,sigma=gaus_estimate(readstar[ic]);
            cent[0],cent[1],sigma1=centriod_psf(readstar[ic]-mean);print(ic,cent[0],cent[1])
            plt.imshow(readstar[ic]**0.1);plt.show()
            readspos[k]=readspos[ic]
            readstar[k]=readstar[ic]
            shiftxy[k]=shiftxy[ic]
            readwgt[k]=readwgt[ic]
            k+=1
    Nobj=k;print("clean Nobj=",Nobj)
    k=0
    for ic in range(Nobj):
        if shiftxy[ic][0]<2 and shiftxy[ic][1]<2:
            k+=1
    spos=np.zeros((k,2),dtype='float')
    instar=np.zeros((k,Ng0,Ng0),dtype='float')
    wgt=np.zeros(instar.shape,dtype='float');print(instar.shape,wgt.shape)
    k=0
    for ic in range(Nobj):
        if shiftxy[ic][0]<3 and shiftxy[ic][1]<3:
            #mean,sigma=gaus_estimate(readstar[ic]);
            #cent[0],cent[1],sigma1=centriod_psf(readstar[ic]-mean);print(ic,cent[0],cent[1])
            spos[k]=readspos[ic]
            instar[k]=readstar[ic]
            wgt[k]=readwgt[ic]
            k+=1
    Nobj=k;print("pos clean Nobj=",Nobj)
    if Nobj>6:
        imodel=np.zeros((Nobj,(Ng0+Nim)*osam,(Ng0+Nim)*osam),dtype='float')
        if Instr=='NIRCAM':
            nrc = webbpsf.NIRCam()
            nrc.filter =Filt
            nrc.detector=detector[Detec]
            nrc.load_wss_opd_by_date(OPD)
        else:
            print("only NIRCAM can be used at present code!")
            sys.exit()
        for ic in range(Nobj):
            print(slicename+" mixpsf imodel",ic+1)
            nrc.detector_position=(spos[ic][0],spos[ic][1])
            psf = nrc.calc_psf(fov_pixels=Ng0+Nim,oversample=osam)
            imodel[ic]=(psf[2].data/np.sum(psf[2].data)).T
        #core function for HybPSF
        PCs,pos,coeff,getNobj,slecstar,webcoeff,slecmask,slecmodel=web_psf_fit(
            instar,
            spos,
            npc,
            gain=float(gains[Detec]),
            method=1,
            wgt=wgt,
            imodel=imodel,
            osam=osam,
            SNRs=10.)
        write_mult_fits(wdirs+slicename+"_PC.fits",PCs,pos,coeff,webcoeff)
        wgt=slecmask;imdoel=slecmodel
        hdu=fits.open(wdirs+slicename+"_PC.fits")
        hdu[0].header['INSTRUME']=Instr
        hdu[0].header['FILTER']=Filt
        hdu[0].header['DETECTOR']=detector[Detec]
        hdu[0].header['DATE']=OPD
        hdu[0].header['OSAM']=osam
        hdu.writeto(wdirs+slicename+"_PC.fits",overwrite=True)
        hdu.close()
        print("The PC data are now stored at path: "+wdirs+"with name: "+slicename+"_PC.fits")
        polyorder=2
        if diagnose==True and getNobj>6:
            if getNobj>10:polyorder=3
            if getNobj>15:polyorder=4
            if getNobj>21:polyorder=5
            rPSF=coeff2psf(pos,coeff,PCs,pos,degrees=polyorder);print("rPSF",rPSF.shape)
            icoeff=polyfit2d(pos,webcoeff[:,0],webcoeff[:,1],pos,3)
            eNg=rPSF.shape[1]
            Ng=int(eNg/osam)
            sNg=60;dpix=int((slecstar.shape[1]-sNg)/2);cut=int(((Ng0+Nim)*osam-eNg)/2)
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
                mean,sigma=gaus_estimate(slecstar[ic])
                instar[ic]-=mean
                snrs=S2N(slecstar[ic],gain=float(gains[Detec]))
                tmp=wgtmap(slecstar[ic],gain=float(gains[Detec]))
                tmp_img=slecstar[ic][dpix:dpix+sNg,dpix:dpix+sNg]
                tmp_wgt=wgt[ic][dpix:dpix+sNg,dpix:dpix+sNg]#mask
                tmp_img*=tmp_wgt
                sums=np.sum(tmp_img);tmp_img/=sums #normalise
                error_map=tmp[dpix:dpix+sNg,dpix:dpix+sNg];error_map/=sums  #weight map
                tmp_wmap=1./(error_map*error_map)*tmp_wgt;
                cent[0],cent[1],sigma=centriod_psf(tmp_img);#print("cx=%f,cy=%f"%(cent[0],cent[1]))
                dx=cent[0]-(sNg/2.+0.5)+1;dy=cent[1]-(sNg/2.+0.5)+1
                tmp_mod=interp_cubic(imodel[ic],sNg,dy,dx,osam)
                tmp_mod*=tmp_wgt;
                mix_mod[ic]*=0
                mix_mod[ic]=rPSF[ic]+imodel[ic,cut:cut+eNg,cut:cut+eNg]*icoeff[ic]
                for i in range(eNg):
                    for j in range(eNg):
                        mix_mod[ic][i][j]=0
                        mix_mod[ic][i][j]=(rPSF[ic][i][j]+imodel[ic][i+cut][j+cut]*icoeff[ic])#mixture model
                        if mix_mod[ic][i][j]<0 : mix_mod[ic][i][j]=0
                for i in range(eNg):
                    for j in range(eNg):
                        if mix_mod[ic][i][j]==0 and i-1>0 and i+1<eNg and j-1>0 and j+1<eNg:
                            mix_mod[ic][i][j]=(mix_mod[ic][i-1][j-1]+mix_mod[ic][i-1][j]+mix_mod[ic][i-1][j+1]+
                                               mix_mod[ic][i][j-1]+mix_mod[ic][i][j+1]+
                                               mix_mod[ic][i+1][j-1]+mix_mod[ic][i+1][j]+mix_mod[ic][i+1][j+1])/8.
                tmp_mix_mod=interp_cubic(mix_mod[ic],sNg,dy,dx,osam)
                mmodel[ic]=tmp_mix_mod
                tmp_mix_mod*=tmp_wgt;
                tmp_img/=np.sum(tmp_img);tmp_mod/=np.sum(tmp_mod);tmp_mix_mod/=np.sum(tmp_mix_mod)
                residu1[ic]=(tmp_img-tmp_mod).T
                residu2[ic]=(tmp_img-tmp_mix_mod).T
                oe1,oe2,or2=size(tmp_img,cent,4.);
                we1,we2,wr2=size(tmp_mod,cent,4.);
                me1,me2,mr2=size(tmp_mix_mod,cent,4.);
                chi1[ic]=residu1[ic]*residu1[ic]*tmp_wmap
                chi2[ic]=residu2[ic]*residu2[ic]*tmp_wmap
                wmodel[ic]=tmp_mod;
                fpw.writelines(str(we1)+'\t'+str(we2)+'\t'+str(wr2)+'\t'+str(spos[ic][0])+'\t'+str(spos[ic][1])+'\t'+str(snrs)+'\n')
                fp.writelines(str(oe1)+'\t'+str(oe2)+'\t'+str(or2)+'\t'+str(spos[ic][0])+'\t'+str(spos[ic][1])+'\t'+str(snrs)+'\n')
                fpm.writelines(str(me1)+'\t'+str(me2)+'\t'+str(mr2)+'\t'+str(spos[ic][0])+'\t'+str(spos[ic][1])+'\t'+str(snrs)+'\n')
            fp.close();fpw.close();fpm.close()
            write_fits(wdirs+slicename+"_webchi.fits",chi1)
            write_fits(wdirs+slicename+"_mixchi.fits",chi2)
            write_fits(wdirs+slicename+"_webresidu.fits",residu1)
            write_fits(wdirs+slicename+"_mixresidu.fits",residu2)
            write_fits(wdirs+slicename+"_wmodel.fits",wmodel)
            write_fits(wdirs+slicename+"_mmodel.fits",mmodel)
            print("the diagnose file are write to "+wdirs)






