import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import ctypes
import webbpsf
import os

from psf_fit import star_shape,size,write_fits,get_filename,\
web_psf_fit,web_psf_rec,int2cent,model2pos,wgtmap,gaus_estimate\
,write_mult_fits,polyfit2d
from C_tools_lib import centriod_psf
from mcmcs import coeff2psf,interp_cubic



nrca_short_detectors = ['NRCA1', 'NRCA2', 'NRCA3', 'NRCA4', 'NRCB1', 'NRCB2', 'NRCB3', 'NRCB4']
nrca_long_detectors = ['NRCA5', 'NRCB5']
F090W_fname="jw02736001001_02101_"
F150W_fname="jw02736001001_02103_"
F200W_fname="jw02736001001_02105_"
F277W_fname="jw02736001001_02101_"
F356W_fname="jw02736001001_02103_"
F444W_fname="jw02736001001_02105_"
sdetect=["nrca1","nrca2","nrca3","nrca4","nrcb1","nrcb2","nrcb3","nrcb4"]
ldetect=["nrcalong","nrcblong"]
ldetect=["nrcalong"]
detector=dict({"nrca1":"nrca1","nrca2":"nrca2","nrca3":"nrca3","nrca4":"nrca4",
               "nrcb1":"nrcb1","nrcb2":"nrcb2","nrcb3":"nrcb3","nrcb4":"nrcb4",
               "nrcalong":"nrca5","nrcblong":"nrcb5"})
filte_name=dict({"F090W":"jw02736001001_02101_","F150W":"jw02736001001_02103_",
                "F200W":"jw02736001001_02105_","F277W":"jw02736001001_02101_",
                "F356W":"jw02736001001_02103_","F444W":"jw02736001001_02105_"})
f2detpath=dict({"F090W":sdetect,"F150W":sdetect,"F200W":sdetect,
                "F277W":ldetect,"F356W":ldetect,"F444W":ldetect})
gains=dict({"nrca1":"2.08","nrca2":"2.02","nrca3":"2.17","nrca4":"2.02",
            "nrcb1":"2.01","nrcb2":"2.14","nrcb3":"1.94","nrcb4":"2.03",
            "nrcalong":"1.84","nrcblong":"1.80"})




osam=2
npc=4

Cam="NIRCam/"
filte="F277W"
dirs="/home/lnie/code/jwstmc/data/"
PCdirs="/home/lnie/code/jwstmc/data/"
TNg=80
ikk=0
Nim=5
for det in f2detpath[filte]:
    gain=float(gains[det])
    stamp_dirs=dirs
    for Tim in range(1):
        get_name=get_filename(dirs,"_00001_nrcalong_cal_pv181_rev4_star.fits")
        ofname=dirs+"jw02736001001_02101_00001_nrcalong_cal.fits"
        hdu=fits.open(ofname);OPD=hdu[0].header['DATE-BEG'];
        print("OPD:",OPD)
        PHOTUJA2=hdu['SCI'].header['PHOTUJA2'];XPOSURE=hdu['SCI'].header['XPOSURE']
        PHOTMJSR=hdu['SCI'].header['PHOTMJSR']
        PIXAR_A2=hdu['SCI'].header['PIXAR_A2']
        tmpname=get_name[0].split("/");l0=len(tmpname);l1=len(tmpname[l0-1]);slicename=tmpname[l0-1][0:l1-len("_star.fits")]
        stamp_name=(stamp_dirs+slicename+"_star.fits")
        cat_name=dirs+slicename+"_star.cat"
        mask_name=dirs+slicename+"_mask.fits"
        print(cat_name,mask_name,filte,det)
        try: 
            cat=np.loadtxt(cat_name)
            hdu=fits.open(stamp_name)
            dpix=int((hdu[0].data.shape[1]-TNg)/2)
            instar=hdu[0].data[:,dpix:TNg+dpix,dpix:TNg+dpix]*XPOSURE*PHOTMJSR;
            print("instar shape:",instar.shape)
            Nobj=instar.shape[0]  
            hdu=fits.open(mask_name)
            inmask=hdu[0].data[:,dpix:TNg+dpix,dpix:TNg+dpix];
        except :
            Nobj=0
            pass
        if Nobj<=6:
            pass
        else:
            #use mixture model......
            #Nobj=int(Nobj*0.8)
            Ng0=instar.shape[1]
            Nh=instar.shape[1]/2
            for ic in range(Nobj):
                instar[ic]=instar[ic].T
                inmask[ic]=inmask[ic].T
            spos=np.zeros((Nobj,2),dtype='float')
            wgt=np.zeros(instar.shape,dtype='float')
            #hdu=fits.open(wgt_name);dwgt=hdu[0].data;dwgt=dwgt.T #read mask maps
            for ic in range(Nobj):
                spos[ic][0]=cat[ic][1];spos[ic][1]=cat[ic][2]
                rx=cat[ic][1]+0.5;ry=cat[ic][2]+0.5
                xl=int(rx-Nh);xh=int(rx+Nh);yl=int(ry-Nh);yh=int(ry+Nh)
                wgt[ic]=inmask[ic][:,:]
                '''if ic <20:
                    plt.subplot(1,3,1);plt.imshow(wgt[ic])
                    plt.subplot(1,3,2);plt.imshow(inmask[ic])
                    plt.subplot(1,3,3);plt.imshow((np.fabs(instar[ic]))**(0.1))
                    plt.show()'''
            imodel=np.zeros((Nobj,(Ng0+Nim)*osam,(Ng0+Nim)*osam),dtype='float')
            nrc = webbpsf.NIRCam()
            nrc.filter =filte
            tmp=det.upper();
            nrc.detector=detector[det].upper()
            nrc.load_wss_opd_by_date(OPD)
            print(filte,detector[det].upper(),OPD)
            for ic in range(Nobj):
                print(slicename+" mixpsf imodel",ic+1,spos[ic][0],spos[ic][1])
                nrc.detector_position=(spos[ic][0],spos[ic][1])
                psf = nrc.calc_psf(fov_pixels=Ng0+Nim,oversample=2)
                imodel[ic]=(psf[0].data/np.sum((psf[0].data))).T
                #imodel[ic]/=np.sqrt(np.sum(imodel[ic]**2))
                #imodel[ic]/=np.sqrt(np.sum(imodel[ic]))
            #wgt=np.ones(instar.shape,dtype='float')
            write_fits("fits/pmodel.fits",imodel)
            PCs,pos,coeff,getNobj,slecstar,webcoeff=web_psf_fit(
                instar,
                spos,
                npc,
                gain,
                method=1,
                wgt=wgt,
                imodel=imodel,
                osam=2,
                SNRs=0.0000001)
            write_mult_fits(PCdirs+slicename+"_PC.fits",PCs,pos,coeff)
            write_fits("fits/PCs.fits",PCs)
            print(PCdirs+slicename+"_PC.fits")
            hdu=fits.open(PCdirs+slicename+"_PC.fits")
            hdu[0].header['filter']=filte
            hdu[0].header['detector']=detector[det].upper()
            hdu[0].header['data']=OPD
            hdu.writeto(PCdirs+slicename+"_PC.fits",overwrite=True)
            hdu.close()
            polyorder=2
            if getNobj>10:polyorder=3
            #if getNobj>15:polyorder=4
            #if getNobj>21:polyorder=5
            print(coeff)
            rPSF=coeff2psf(pos,coeff,PCs,pos,polyorder);print("1",rPSF.shape)
            icoeff=polyfit2d(pos,webcoeff[:,0],webcoeff[:,1],pos,polyorder);print("2")
            #icoeff=webcoeff[:,0]#polyfit2d(pos,webcoeff[:,0],webcoeff[:,1],pos,polyorder=1)
            eNg=rPSF.shape[1]
            Ng=int(eNg/osam)
            '''for ic in range(getNobj):
                for i in range(eNg):
                    for j in range(eNg):
                        if np.fabs(rPSF[ic][i][j])<10**(-4):rPSF[ic][i][j]=0'''
            sNg=60;dpix=int((slecstar.shape[1]-sNg)/2);cut=int(((Ng0+Nim)*osam-eNg)/2)
            cent=np.zeros(2,dtype='float')
            mix_mod=np.zeros(rPSF.shape,dtype='float')
            rrPSF=np.zeros((getNobj,sNg,sNg),dtype='float')
            chi1=np.zeros((getNobj,sNg,sNg))
            chi2=np.zeros((getNobj,sNg,sNg))
            residu1=np.zeros((getNobj,sNg,sNg))
            residu2=np.zeros((getNobj,sNg,sNg))
            fp=open(stamp_dirs+slicename+"_oshape.cat","w")
            fpw=open(stamp_dirs+slicename+"_wshape.cat","w")
            fpm=open(stamp_dirs+slicename+"_mshape.cat","w")
            for ic in range(getNobj):
                mean,sigma=gaus_estimate(slecstar[ic]);
                print(ic,mean,sigma)
                #mean,sigma=gaus_estimate(rPSF[ic]);rPSF[ic]-=mean
                slecstar[ic]-=mean
                #slecstar[ic]*=wgt[ic]
                slecstar[ic]/=np.sum(slecstar[ic])
                #snrs=S2N(instar[ic],gain)
                tmp=wgtmap(slecstar[ic])
                tmp_img=slecstar[ic][dpix:dpix+sNg,dpix:dpix+sNg]
                tmp_wgt=wgt[ic][dpix:dpix+sNg,dpix:dpix+sNg]#mask
                tmp_img*=tmp_wgt
                #tmp_img/=np.sum(tmp_img) #normalise
                error_map=tmp[dpix:dpix+sNg,dpix:dpix+sNg];#error_map/=sums  #weight map
                tmp_wmap=1./(error_map*error_map)*tmp_wgt;
                cent[0],cent[1],sigma=centriod_psf(tmp_img);
                print(cent[0],cent[1])
                dx=cent[0]-(sNg/2.+0.5)+1;dy=cent[1]-(sNg/2.+0.5)+1
                rrPSF[ic]=interp_cubic(rPSF[ic],sNg,dy,dx,osam=2)
                tmp_mod=interp_cubic(imodel[ic],sNg,dy,dx,osam=2)
                #print("sum:", np.sum(rrPSF[ic]))
                #mean,sigma=gaus_estimate(rrPSF[ic]);rrPSF[ic]-=mean
                rrPSF[ic]=rrPSF[ic].T;
                tmp_mod*=tmp_wgt;
                frac=(1.-np.sum(imodel[ic,cut:cut+eNg,cut:cut+eNg]))
                #rPSF=rPSF/np.sum(rPSF)*frac
                #tmp_mod/=np.sum(tmp_mod)
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
                tmp_mix_mod=interp_cubic(mix_mod[ic],sNg,dy,dx,osam)
                tmp_mix_mod*=tmp_wgt;
                print(np.sum(tmp_img),np.sum(tmp_mix_mod))
                #mean,sigma=gaus_estimate(tmp_mix_mod);tmp_mix_mod-=mean
                tmp_img/=np.sum(tmp_img);tmp_mod/=np.sum(tmp_mod);tmp_mix_mod/=np.sum(tmp_mix_mod)
                residu1[ic]=(tmp_img-tmp_mod).T
                residu2[ic]=(tmp_img-tmp_mix_mod).T
                oe1,oe2,or2=size(tmp_img,cent,4.);
                we1,we2,wr2=size(tmp_mod,cent,4.);
                me1,me2,mr2=size(tmp_mix_mod,cent,4.);
                fpw.writelines(str(we1)+'\t'+str(we2)+'\t'+str(wr2)+'\t'+str(spos[ic][0])+'\t'+str(spos[ic][1])+'\n')
                fp.writelines(str(oe1)+'\t'+str(oe2)+'\t'+str(or2)+'\t'+str(spos[ic][0])+'\t'+str(spos[ic][1])+'\n')
                fpm.writelines(str(me1)+'\t'+str(me2)+'\t'+str(mr2)+'\t'+str(spos[ic][0])+'\t'+str(spos[ic][1])+'\n')
                chi1[ic]=residu1[ic]*residu1[ic]*tmp_wmap
                chi2[ic]=residu2[ic]*residu2[ic]*tmp_wmap
            write_fits("fits/rPSF.fits",rrPSF)
            fp.close();fpw.close();fpm.close()
            write_fits(stamp_dirs+slicename+"_webchi.fits",chi1)
            write_fits(stamp_dirs+slicename+"_mixchi.fits",chi2)
            write_fits(stamp_dirs+slicename+"_webresidu.fits",residu1)
            write_fits(stamp_dirs+slicename+"_mixresidu.fits",residu2)  
            write_fits("fits/ostars.fits",instar)
            oshape=np.loadtxt(stamp_dirs+slicename+"_oshape.cat")
            wshape=np.loadtxt(stamp_dirs+slicename+"_mshape.cat")
            de1=oshape[:,0]-wshape[:,0];de2=oshape[:,1]-wshape[:,1];dr=(oshape[:,2]**0.5-wshape[:,2]**0.5)
            print(np.mean(oshape[:,0]),np.std(oshape[:,0]),np.mean(oshape[:,1]),np.std(oshape[:,1]),np.mean(oshape[:,2]),np.std(oshape[:,2]))
            print(np.mean(de1),np.std(de1),np.mean(de2),np.std(de2),np.median(dr),np.std(dr),np.median((oshape[:,2]-wshape[:,2])/oshape[:,2]))










        

