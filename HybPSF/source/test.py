import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import ctypes
import webbpsf
import os
from scipy.stats import sigmaclip
from astropy.stats import SigmaClip

from psf_fit import star_shape,write_fits,get_filename,\
web_psf_fit,web_psf_rec,int2cent,model2pos,wgtmap,gaus_estimate\
,write_mult_fits,polyfit2d
from C_tools_lib import centriod_psf
from mcmcs import coeff2psf,interp_cubic

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



#rd=0.9time
osam=8
npc=10
#npc=12

Cam="NIRCam/"
filte="F444W"
odirs="/data/SMACS0723_NIRCam_pipeline1.8.1/"+filte+'/'
dirs="/home/lnie/data/JWST_rev/"+Cam+filte+'/'
PCdirs="/home/lnie/data/PC/JWST/"+Cam+filte+'/'
TNg=80
ikk=0
Nim=5
cent=np.zeros(2,dtype='float')
for det in f2detpath[filte]:
    gain=float(gains[det])
    stamp_dirs=dirs+det+'/'+"catalogue/star_stamps/"
    for Tim in range(1):
        expos=Tim+1
        indx=str(int(expos/10000))+str(int(expos/1000)%10)+str(int(expos/100)%10)+str(int(expos/10)%10)+str(expos%10)
        ofname=(odirs+filte_name[filte]+indx+'_'+det+"_cal.fits");print(ofname)
        stamp_name=(stamp_dirs+filte_name[filte]+indx+'_'+det+"_cal_pv181_rev4_star.fits")
        tmpname=stamp_name.split("/");l0=len(tmpname);l1=len(tmpname[l0-1]);slicename=tmpname[l0-1][0:l1-len("_star.fits")]
        hdu=fits.open(ofname);OPD=hdu[0].header['DATE-BEG'];
        print("OPD:",OPD)
        PHOTUJA2=hdu['SCI'].header['PHOTUJA2'];XPOSURE=hdu['SCI'].header['XPOSURE']
        PHOTMJSR=hdu['SCI'].header['PHOTMJSR']
        PIXAR_A2=hdu['SCI'].header['PIXAR_A2']
        stamp_name=(stamp_dirs+slicename+"_star.fits")
        cat_name=dirs+det+'/'+"catalogue/star/"+slicename+"_star.cat"
        mask_name=dirs+det+'/'+"catalogue/star_stamps/"+slicename+"_mask.fits"
        print(cat_name,mask_name,filte,det)
        try: 
            cat=np.loadtxt(cat_name)
            hdu=fits.open(stamp_name)
            dpix=int((hdu[0].data.shape[1]-TNg)/2)
            readstar=hdu[0].data[:,dpix:TNg+dpix,dpix:TNg+dpix]*XPOSURE*PHOTMJSR;
            print("instar shape:",readstar.shape)
            Nobj=readstar.shape[0]  
            hdu=fits.open(mask_name)
            readmask=hdu[0].data[:,dpix:TNg+dpix,dpix:TNg+dpix];
        except :
            Nobj=0
            pass
        if Nobj<7:
            pass
        else:
            #use mixture model......
            #Nobj=int(Nobj*0.8)
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
                    #plt.imshow(readstar[ic]**0.1);plt.show()
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
                if shiftxy[ic][0]<2 and shiftxy[ic][1]<2:
                    #mean,sigma=gaus_estimate(readstar[ic]);
                    #cent[0],cent[1],sigma1=centriod_psf(readstar[ic]-mean);print(ic,cent[0],cent[1])
                    spos[k]=readspos[ic]
                    instar[k]=readstar[ic]
                    wgt[k]=readwgt[ic]
                    k+=1
            Nobj=k;print("pos clean Nobj=",Nobj)
            if Nobj>6:
                try:
                    imodel=fits.open(dirs+det+'/'+"catalogue/star_stamps/"+slicename+"_imodel3_latest.fits")[0].data
                except:
                    print("generate wmodel",Nobj)
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
                        psf = nrc.calc_psf(fov_pixels=Ng0+Nim,oversample=osam)
                        imodel[ic]=(psf[2].data/np.sum((psf[2].data))).T
                    write_fits(dirs+det+'/'+"catalogue/star_stamps/"+slicename+"_imodel_latest.fits",imodel)
                PCs,pos,coeff,getNobj,slecstar,webcoeff,slecwgt,slecimodel=web_psf_fit(
                    instar,
                    spos,
                    npc,
                    gain,
                    method=1,
                    wgt=wgt,
                    imodel=imodel,
                    osam=osam,
                    SNRs=10.)
                write_mult_fits(PCdirs+slicename+"_PC_latest.fits",PCs,pos,coeff,webcoeff)
                wgt=slecwgt;imodel=slecimodel
                #write_fits("fits/PCs.fits",PCs)
                print(PCdirs+slicename+"_PC_latest.fits")
                hdu=fits.open(PCdirs+slicename+"_PC_latest.fits")
                hdu[0].header['FILTER']=filte
                hdu[0].header['DETECTOR']=detector[det].upper()
                hdu[0].header['DATE']=OPD
                hdu.writeto(PCdirs+slicename+"_PC_latest.fits",overwrite=True)
                hdu.close()
                polyorder=2
                if getNobj>10:polyorder=3
                if getNobj>15:polyorder=4
                if getNobj>21:polyorder=5
                if getNobj>5 :
                    rPSF=coeff2psf(pos,coeff,PCs,pos,degrees=polyorder);print("rPSF",rPSF.shape)
                    icoeff=polyfit2d(pos,webcoeff[:,0],webcoeff[:,1],pos,3)
                    eNg=rPSF.shape[1]
                    Ng=int(eNg/osam)
                    sNg=60;dpix=int((slecstar.shape[1]-sNg)/2);cut=int(((Ng0+Nim)*osam-eNg)/2)
                    mix_mod=np.zeros(rPSF.shape,dtype='float')
                    rrPSF=np.zeros((getNobj,sNg,sNg),dtype='float')
                    chi1=np.zeros((getNobj,sNg,sNg))
                    chi2=np.zeros((getNobj,sNg,sNg))
                    residu1=np.zeros((getNobj,sNg,sNg))
                    residu2=np.zeros((getNobj,sNg,sNg))
                    fp=open(stamp_dirs+slicename+"_mcoshape1_latest.cat","w")
                    fpw=open(stamp_dirs+slicename+"_mcwshape1_latest.cat","w")
                    fpm=open(stamp_dirs+slicename+"_mcmshape1_latest.cat","w")
                    nrc = webbpsf.NIRCam()
                    nrc.filter =filte
                    tmp=det.upper();
                    nrc.detector=detector[det].upper()
                    nrc.load_wss_opd_by_date(OPD)
                    for ic in range(getNobj):
                        mean,sigma=gaus_estimate(slecstar[ic]);
                        #plt.imshow((np.fabs(slecstar[ic]))**0.5)
                        print(ic,mean,sigma)
                        #mean,sigma=gaus_estimate(rPSF[ic]);rPSF[ic]-=mean
                        slecstar[ic]-=mean
                        snrs=S2N(slecstar[ic],gain)
                        #slecstar[ic]*=wgt[ic]
                        tmp=wgtmap(slecstar[ic])
                        slecstar[ic]/=np.sum(slecstar[ic])
                        tmp_img=slecstar[ic][dpix:dpix+sNg,dpix:dpix+sNg]
                        tmp_wgt=wgt[ic][dpix:dpix+sNg,dpix:dpix+sNg]#mask
                        tmp_img*=tmp_wgt
                        #tmp_img/=np.sum(tmp_img) #normalise
                        error_map=tmp[dpix:dpix+sNg,dpix:dpix+sNg];#error_map/=sums  #weight map
                        tmp_wmap=1./(error_map*error_map)*tmp_wgt;
                        cent[0],cent[1],sigma=centriod_psf(tmp_img);
                        print(cent[0],cent[1])
                        dx=cent[0]-(sNg/2.+0.5)+1;dy=cent[1]-(sNg/2.+0.5)+1
                        try:
                            rrPSF[ic]=interp_cubic(rPSF[ic],sNg,dy,dx,osam=osam)
                            #tmp_mod=interp_cubic(imodel[ic],sNg,dy,dx,osam=osam)
                            dx,dy=cent[0]-sNg/2+0.5,cent[1]-sNg/2+0.5
                            nrc.detector_position=(spos[ic][0],spos[ic][1])
                            nrc.options['source_offset_x'] = (dx*nrc.pixelscale)
                            nrc.options['source_offset_y'] = (dy*nrc.pixelscale)
                            psf = nrc.calc_psf(fov_pixels=sNg)
                            tmp_mod=(psf[3].data/np.sum((psf[3].data))).T
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
                                    if mix_mod[ic][i][j]==0 and i-1>0 and i+1<eNg and j-1>0 and j+1<eNg:
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
                            print(oe1,oe2,or2)
                            print(me1,me2,mr2)
                            fpw.writelines(str(we1)+'\t'+str(we2)+'\t'+str(wr2)+'\t'+str(spos[ic][0])+'\t'+str(spos[ic][1])+'\t'+str(snrs)+'\n')
                            fp.writelines(str(oe1)+'\t'+str(oe2)+'\t'+str(or2)+'\t'+str(spos[ic][0])+'\t'+str(spos[ic][1])+'\t'+str(snrs)+'\n')
                            fpm.writelines(str(me1)+'\t'+str(me2)+'\t'+str(mr2)+'\t'+str(spos[ic][0])+'\t'+str(spos[ic][1])+'\t'+str(snrs)+'\n')
                            chi1[ic]=residu1[ic]*residu1[ic]*tmp_wmap
                            chi2[ic]=residu2[ic]*residu2[ic]*tmp_wmap
                        except:
                            pass
                    #write_fits("fits/rPSF.fits",rrPSF)
                    fp.close();fpw.close();fpm.close()
                    write_fits(stamp_dirs+slicename+"_webchi_latest.fits",chi1)
                    write_fits(stamp_dirs+slicename+"_mixchi_latest.fits",chi2)
                    write_fits(stamp_dirs+slicename+"_webresidu_latest.fits",residu1)
                    write_fits(stamp_dirs+slicename+"_mixresidu_latest.fits",residu2)  
                    #write_fits("fits/ostars.fits",instar)
                    oshape=np.loadtxt(stamp_dirs+slicename+"_mcoshape1_latest.cat")
                    wshape=np.loadtxt(stamp_dirs+slicename+"_mcmshape1_latest.cat")
                    de1=oshape[:,0]-wshape[:,0];de2=oshape[:,1]-wshape[:,1];dr=(oshape[:,2]**0.5-wshape[:,2]**0.5)
                    print(np.mean(oshape[:,0]),np.mean(oshape[:,1]),np.mean(oshape[:,2]))
                    print(np.mean(de1),np.std(de1),np.mean(de2),np.std(de2),np.median(dr),np.std(dr),np.mean((oshape[:,2]-wshape[:,2])/oshape[:,2]))
                else:
                    pass
            else:
                pass

        


