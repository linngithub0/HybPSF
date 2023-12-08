import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os
import ctypes
from scipy import linalg,interpolate
import sys
import platform
from scipy.stats import sigmaclip
import webbpsf
#out

def _get_filename(file_dir, str1_match="default_string", str2_match="default_string"):
    fname = []
    if str2_match == "default_string":
        for list in os.listdir(file_dir):
            if os.path.isfile(os.path.join(file_dir, list)):
                if str(str1_match) in list:
                    fname.append(list)
    else:
        for list in os.listdir(file_dir):
            if os.path.isfile(os.path.join(file_dir, list)):
                if str(str1_match) in list and str(str2_match) in list:
                    fname.append(list)

    return fname

libdir = os.path.dirname(os.path.abspath(__file__))
tmp = os.path.join(libdir, "source")
fnlib = _get_filename(tmp, '.so', "libpsffit")[0]
psflibdir = os.path.join(libdir, "source", fnlib)
_libpsffit = ctypes.CDLL(psflibdir)
tmp = os.path.join(libdir, "staridf")
fnlib = _get_filename(tmp, ".so", "star_extrac")[0]
starlibdir = os.path.join(libdir, "staridf", fnlib)
_star_extrac = ctypes.CDLL(starlibdir)
fpoint = ctypes.POINTER(ctypes.c_float)
dpoint = ctypes.POINTER(ctypes.c_double)




fpoint=ctypes.POINTER(ctypes.c_float)
dpoint=ctypes.POINTER(ctypes.c_double)





detector=dict({"NRCA1":"NRCA1","NRCA2":"NRCA2","NRCA3":"NRCA3","NRCA4":"NRCA3",
               "NRCB1":"NRCB1","NRCB2":"NRCB2","NRCB3":"NRCB3","NRCB4":"NRCB4",
               "NRCALONG":"NRCA5","NRCABLONG":"NRCB5"})
gains=dict({"NRCA1":"2.08","NRCA2":"2.02","NRCA3":"2.17","NRCA4":"2.02",
            "NRCB1":"2.01","NRCB2":"2.14","NRCB3":"1.94","NRCB4":"2.03",
            "NRCALONG":"1.84","NRCABLONG":"1.80"})
#the gain information are obtained from Table 2 of link:
#https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-instrumentation/nircam-detector-overview/nircam-detector-performance


cat_header=["#   1 NUMBER                 Running object number\n",                                  
"#   2 X_IMAGE                Object position along x                                      [pixel]\n",
"#   3 Y_IMAGE                Object position along y                                      [pixel]\n",
"#   4 FWHM_IMAGE             FWHM assuming a gaussian core                                [pixel]\n",
"#   5 FLUX_RADIUS            Fraction-of-light radii                                      [pixel]\n",
"#   6 FLUX_AUTO              Flux within a Kron-like elliptical aperture                  [count]\n",
"#   7 FLUXERR_AUTO           RMS error for AUTO flux                                      [count]\n",
"#   8 FLAGS                  Extraction flags                                                    \n"                                          
"#   9 CLASS_STAR             S/G classifier output\n"]

def _Image2PSFmodel(fullfitsname,wdirs,NPCs,SNR_thresh=20,osam=1,
    stamp_size=80,magleft=-15.,magright=-12.,big=0):
    """
    select the PSF star candidates from a given fitsfile file

    Parameters:
    ----------
    fullfitsname : str
        The full name of the fits file, defaulting to 'fitsfile/name'.
    wdirs: str
        the path to save the selected star stamps
    NPCs: int
        number of PCs, larger than 4 are reommended
    SNR_thresh: float
        signal to noise ratio threshold for the star images,which will used in PCs calculation
    osam : int
        Oversampling factor of PSF model, >1 is recommended.
    stamp_size: int
        the pxiel size of the cutout star stamp images:stamp_size*stamp_size
    magleft: float
        a threshold to select the star images,which are defined log10(flux)
    magright: float
        a threshold to select the star images,which are defined log10(flux)
    big: int
        a chose for PSF reconstruction situation
        0 for samll source
        1 for very bright source

    Returns:
    ----------
    str: 
        the file name of the PSF moddel file.
    """
    if sys.platform.startswith('win'):
        tmpname=fullfitsname.split('\\');
    else:
        tmpname=fullfitsname.split('/');
    l0=len(tmpname);
    fdirs='/'
    for i in range(l0-1):
        fdirs=os.path.join(fdirs,tmpname[i])
    print(fdirs)
    fitsname=tmpname[l0-1]
    try:
        print("_pick_psfstars")
        stamp_name=_pick_psfstars(fdirs,fitsname,wdirs,stamp_size,magleft=magleft,magright=magright)
        stars=fits.open(stamp_name)[0].data
        if stars.shape[0]>NPCs:
            PCname=_get_PC(wdirs,stamp_name,NPCs=NPCs,method=1,osam=osam,SNR_thresh=SNR_thresh,big=big)
            return(PCname)
        else:
            return(0)
    except:
        return(0)



def _pick_psfstars(fdirs,fitsname,wdirs,stamp_size=80,magleft=-15.,magright=-12.):
    """
    select the PSF star candidates from a given fitsfile file

    Parameters:
    ----------
    fdirs: str
        the file path of the fitsfile
    fitsname: str
        the fitsfile name of the fits file
    wdirs: str
        the path to save the selected star stamps
    stamp_size: int
        the pxiel size of the cutout star stamp images:stamp_size*stamp_size
    magleft: float
        a threshold to select the star images,which are defined log10(flux)
    magright: float
        a threshold to select the star images,which are defined log10(flux)

    Returns:
    ----------
    stamp_name：str
        the file name of star stamps
    cat_nam: str
        the catalogue name of star images
    """
    if os.system("which sex") == 256 and os.system("which source-extractor") == 256:
        raise OSError('No sex SHELL command found! Please install sextractor')
    fname=os.path.join(fdirs,fitsname)
    if sys.platform.startswith('win'):
        tmpname=fname.split("\\");l0=len(tmpname);l1=len(tmpname[l0-1])
    else:
        tmpname=fname.split("/");l0=len(tmpname);l1=len(tmpname[l0-1])
    slicename=tmpname[l0-1][0:l1-len(".fits")]
    hdu=fits.open(fname)
    data=hdu[1].data
    IMGNAXIS1=hdu[1].header['NAXIS1']
    IMGNAXIS2=hdu[1].header['NAXIS2']
    hdu=fits.open(fname);OPD=hdu[0].header['DATE-BEG'];
    print("OPD:",OPD)
    PHOTUJA2=hdu['SCI'].header['PHOTUJA2'];XPOSURE=hdu['SCI'].header['XPOSURE']
    PHOTMJSR=hdu['SCI'].header['PHOTMJSR']
    PIXAR_A2=hdu['SCI'].header['PIXAR_A2']
    INSTRUME=hdu[0].header['INSTRUME']
    DETECTOR=hdu[0].header['DETECTOR']
    FILTER=hdu[0].header['FILTER']
    gain=gains[DETECTOR]
    #find the catalogue using cold+hot mode()
    coldcat=os.path.join(wdirs,"cold.cat")
    libdir = os.path.dirname(os.path.abspath(__file__))
    libsex=os.path.join(libdir,"config","default.sex")
    libconv=os.path.join(libdir,"config","default.conv")
    libnnw=os.path.join(libdir,"config","default.nnw")
    libparam=os.path.join(libdir,"config","default.param")
    print("search source based on cold mode.")
    coldcmd=("sex "+fname+" -c "+libsex+" -CATALOG_NAME "+coldcat+" -DETECT_MINAREA 30 -DETECT_THRESH 2.2 \
        -DEBLEND_NTHRESH 32 -DEBLEND_MINCONT 0.04 -BACK_SIZE 100 -BACK_FILTERSIZE 3 \
        -BACKPHOTO_TYPE LOCAL -FILTER_NAME "+libconv+" -STARNNW_NAME "+libnnw+" -PARAMETERS_NAME "+libparam)
    print("coldcmd",coldcmd)
    os.system(coldcmd)
    if os.system("which source-extractor") != 256:
        coldcmd = ("source-extractor "+fname+" -c "+libsex+" -CATALOG_NAME "+coldcat+" -DETECT_MINAREA 30 -DETECT_THRESH 2.2 \
        -DEBLEND_NTHRESH 32 -DEBLEND_MINCONT 0.04 -BACK_SIZE 100 -BACK_FILTERSIZE 3 \
        -BACKPHOTO_TYPE LOCAL -FILTER_NAME "+libconv+" -STARNNW_NAME "+libnnw+" -PARAMETERS_NAME "+libparam)
        print("coldcmd", coldcmd)
        os.system(coldcmd)
    if os.system("which sex") != 256:
        coldcmd = ("sex "+fname+" -c "+libsex+" -CATALOG_NAME "+coldcat+" -DETECT_MINAREA 30 -DETECT_THRESH 2.2 \
        -DEBLEND_NTHRESH 32 -DEBLEND_MINCONT 0.04 -BACK_SIZE 100 -BACK_FILTERSIZE 3 \
        -BACKPHOTO_TYPE LOCAL -FILTER_NAME "+libconv+" -STARNNW_NAME "+libnnw+" -PARAMETERS_NAME "+libparam)
        print("coldcmd", coldcmd)
        os.system(coldcmd)
    else:
        raise OSError('No sex SHELL command found! Please install sextractor')
    '''
    os.system(coldcmd)
    print("search source based on hot mode.")
    hotcat=os.path.join(wdirs,"hot.cat")
    hotcmd=("sex "+fname+" -c "+libdir+" -CATALOG_NAME "+hotcat+" -DETECT_MINAREA 18 -DETECT_THRESH 1.2 \
        -DEBLEND_NTHRESH 32 -DEBLEND_MINCONT 0.065 -BACK_SIZE 100 -BACK_FILTERSIZE 3 -BACKPHOTO_TYPE LOCAL")
    os.system(hotcmd)
    '''
    try:
        #hot_cat=np.loadtxt(hotcat);print("hot:",hot_cat.shape)
        #cold_cat=hot_cat
        cold_cat=np.loadtxt(coldcat);
        indx=np.where(cold_cat[:,5]>0)
        cold_cat=cold_cat[indx[0]]#select only have flux larger than 0
        print("total slected source number:",cold_cat.shape[0])
        fpn=os.path.join(wdirs,slicename+".cat")
        fp=open(fpn,"w")
        for line in cat_header:
            fp.writelines(line)
        for it in range(cold_cat.shape[0]):
            fp.writelines(str(it+1)+'\t'+str(cold_cat[it][1])+'\t'\
                +str(cold_cat[it][2])+'\t'+str(cold_cat[it][3])+'\t'\
                +str(cold_cat[it][4])+'\t'+str(cold_cat[it][5])+'\t'\
                +str(cold_cat[it][6])+'\t'+str(cold_cat[it][7])+'\t'\
                +str(cold_cat[it][8])+'\n')
        fp.close()
        r50=cold_cat[:,4];fwhm=cold_cat[:,3];mag=-2.5*np.log10(cold_cat[:,5])
        fig=plt.figure()
        plt.plot(mag,r50,'.',c='blue',alpha=0.3)
        plt.plot(mag,fwhm,'.',c='red',alpha=0.5)
        plt.ylim(-2,30)
        fgn=os.path.join(wdirs,slicename+'_cat.pdf')
        #plt.savefig(fgn)
        indx=np.where(cold_cat[:,8]>0.1)
        class_star=cold_cat[indx] #read the candidate star from sextractor star class parameter
        flux=class_star[:,5] #read the flux of class_star
        mag=-2.5*np.log10(flux)
        #bound=(np.max(mag)-np.min(mag))*0.3/2.
        bound1=(np.max(mag)-np.min(mag))*0.3
        bound2=(np.max(mag)-np.min(mag))*0.1
        magleft=np.min(mag)+bound1
        magright=np.max(mag)-bound2
        print("original mag range:",np.min(mag),np.max(mag),"magleft and magright=",magleft,magright)

        print("try find star in : %s"%(wdirs+slicename+".cat"))

        if sys.platform.startswith('win'):
            fdirs=fdirs+"\\"
            wdirs=wdirs+"\\"
        else:
            fdirs=fdirs+"/"
            wdirs=wdirs+'/'
        _sex_star(wdirs,slicename+".cat",fdirs,tmpname[l0-1],wdirs,str(stamp_size),str(magleft),str(magright),str(gain))
        
        #os.system(staridf)
        cat=os.path.join(wdirs,slicename+".cat")
        cold_cat=np.loadtxt(cat)
        cat_name=os.path.join(wdirs,"catalogue","star",slicename+"_star"+".cat")
        #do sigma cliping
        try:
            stamp_name = os.path.join(wdirs, "catalogue",
                                      "star_stamps", "".join([slicename, "_star.fits"]))
            mask_name = os.path.join(wdirs, "catalogue",
                                      "star_stamps", "".join([slicename, "_mask.fits"]))
            rcat = np.loadtxt(cat_name)
            hdu = fits.open(stamp_name)
            stars = hdu[0].data
            hdu = fits.open(mask_name)
            masks = hdu[0].data
            hdu.close()
            # sigmaclip
            Nobj = stars.shape[0]
            print("raw star number is:", Nobj)
            Ng0 = stars.shape[1]
            Nh = int(Ng0/2)
            cent = np.zeros(2, dtype='float')
            shiftxy = []
            shapes = []
            for ic in range(Nobj):
                mean, sigma = _gaus_estimate(stars[ic])
                cent[0], cent[1], sigma1 = _centriod_psf(stars[ic]-mean)
                oe1, oe2, or2 = _size(stars[ic], cent, 5.)
                mage = (oe1**2+oe2**2)**0.5
                shiftxy.append([np.fabs(cent[0]-Ng0/2+0.5), np.fabs(cent[1]-Ng0/2+0.5)])
                shapes.append([mage, or2])
            shapes = np.array(shapes)
            shiftxy = np.array(shiftxy)
            cleane, low, high = sigmaclip(shapes[:, 0], 3., 3.)
            k = 0
            for ic in range(Nobj):
                if (k < cleane.shape[0] and shapes[ic][0] == cleane[k]):
                    rcat[k] = rcat[ic]
                    stars[k] = stars[ic]
                    masks[k] = masks[ic]
                    shiftxy[k] = shiftxy[ic]
                    k += 1
            Nobj = k
            print("clean Nobj=", Nobj)
            k = 0
            for ic in range(Nobj):
                if shiftxy[ic][0] < 1.5 and shiftxy[ic][1] < 1.5:
                    k += 1
            spos = np.zeros((k, rcat.shape[1]), dtype='float')
            instar = np.zeros((k, Ng0, Ng0), dtype='float')
            inmask = np.zeros((k, Ng0, Ng0), dtype='float')
            k = 0
            for ic in range(Nobj):
                if shiftxy[ic][0] < 1.5 and shiftxy[ic][1] < 1.5:
                    spos[k] = rcat[ic]
                    instar[k] = stars[ic]
                    inmask[k] = masks[ic]
                    k += 1
            Nobj = k
            print("pos clean Nobj=", Nobj)
            # plot mag-size
            r50 = spos[:, 4]
            fwhm = spos[:, 3]
            mag = -2.5*np.log10(spos[:, 5])
            plt.plot(mag, r50, '*', c='blue', alpha=0.3, label='r50')
            plt.plot(mag, fwhm, '*', c='red', alpha=0.5, label='fwhm')
            plt.savefig(fgn)
            #hdu=fits.open(fname);
            #rcat=np.loadtxt(cat_name);
            #print(spos.shape,instar.shape,inmask.shape)
            #stamp_name=os.path.join(wdirs,slicename+"_star.fits");print(stamp_name)
            whdu = fits.PrimaryHDU(instar)
            whdu1 = fits.ImageHDU(spos)
            whdu2 = fits.ImageHDU(inmask)
            #whdu3 = fits.ImageHDU(inmask)
            hdul = fits.HDUList([whdu, whdu1, whdu2])
            hdul[0].header['GAIN']=str(gain);
            hdul[0].header['NAXIS1']=IMGNAXIS1
            hdul[0].header['NAXIS2']=IMGNAXIS2
            hdul[0].header['DATE-BEG']=OPD;
            hdul[0].header['PHOTUJA2']=PHOTUJA2;
            hdul[0].header['XPOSURE']=XPOSURE
            hdul[0].header['PHOTMJSR']=PHOTMJSR
            hdul[0].header['PIXAR_A2']=PIXAR_A2
            hdul[0].header['PIXAR_A2']=PIXAR_A2
            hdul[0].header['INSTRUME']=INSTRUME
            hdul[0].header['DETECTOR']=DETECTOR
            hdul[0].header['FILTER']=FILTER
            stamp_name = os.path.join(wdirs, "".join([slicename, "_star.fits"]))
            hdul.writeto(stamp_name,overwrite=True)
            hdul.close()
            return(stamp_name)
        except:
            stamp_name = os.path.join(wdirs, "".join([slicename, "_star.fits"]));
            print(stamp_name+" does not have enough stars")
            pass
    except:
        print("file open failed, please check this file %s"%(fname))
        pass

    

def _get_PC(wdirs,stamp_name,NPCs=10,method=1,osam=1,SNR_thresh=80,big=0):
    '''
    calculate the principal components(PCs) from given star images

    Parameters:
    -------
    stamp_name：str
        the file name of star stamps, which is the return of pick_psfstars
    cat_nam: str
        the catalogue name of star images, which is the return of pick_psfstars
    NPCs: int
        number of PCs, larger than 4 are reommended
    method: string
        reconstruction method used 'iSPCA', 'SPCA', 'EMPCA' or 'PCA'
    S2N: float
        signal to noise ratio threshold for the star images,which will used in PCs calculation
    
    Returns:
    ----------
    PCs: numpy 3d array
        the principal components(PCs)
    spos: numpy 2d array
        the selected star image positions aranged as spos[ic][0]->x,spos[ic][1]->y, 
        this will be used for the PSF field variation modelings
    coeff: numpy 2d array
        the coefficients and the corresponding errors of PCs
    '''
    #mask_name=stamp_name.replace("_star.fits","_mask.fits")
    #print("stamp_name:",stamp_name)
    hdu=fits.open(stamp_name);#print(hdu.info())
    gain=float(hdu[0].header['GAIN'])
    instars=hdu[0].data;
    cat=hdu[1].data
    readmask=hdu[2].data;#print(instars.shape,cat.shape,readmask.shape)
    dpix=5
    psf_size=instars.shape[1]-dpix
    cent=np.zeros(2,dtype='float')
    if sys.platform.startswith('win'):
        tmpname=stamp_name.split("\\");l0=len(tmpname);l1=len(tmpname[l0-1]);slicename=tmpname[l0-1][0:l1-len('.fits')]
    else:
        tmpname=stamp_name.split("/");l0=len(tmpname);l1=len(tmpname[l0-1]);slicename=tmpname[l0-1][0:l1-len('.fits')]
    '''
    fig=plt.figure()
    snrs=[]
    for star in instars:
        snr=_S2N(star,gain=1.1)
        snrs.append(snr)
    plt.plot(snrs,'-')
    plt.savefig(os.path.join(wdirs,slicename+"_snrs.pdf"))
    '''
    #further stamps selection for PSF modeling
    try: 
        #cat=np.loadtxt(cat_name)
        #hdu=fits.open(stamp_name)
        OPD=hdu[0].header['DATE-BEG'];
        print("OPD:",OPD)
        PHOTUJA2=hdu[0].header['PHOTUJA2'];XPOSURE=hdu[0].header['XPOSURE']
        PHOTMJSR=hdu[0].header['PHOTMJSR']
        PIXAR_A2=hdu[0].header['PIXAR_A2']
        FILTER=hdu[0].header['FILTER']
        DETECTOR=hdu[0].header['DETECTOR']
        readstar=hdu[0].data
        #readstar*=(XPOSURE*PHOTMJSR);
        print("instar shape:",readstar.shape)
        Nobj=readstar.shape[0]  
        #hdu=fits.open(mask_name)
        #readmask=hdu[0].data
    except :
        Nobj=0
        print("file: %s not found"%(stamp_name))
        pass
    if Nobj<7:
        pass
        print("star number: %d is not enough to construct model"%(Nobj))
    else:
        #use mixture model......
        Ng0=readstar.shape[1]
        Nh=readstar.shape[1]/2
        Nim=5
        for ic in range(Nobj):
            readstar[ic]=readstar[ic].T
            readmask[ic]=readmask[ic].T
        spos=cat[:,1:3]
        instar=readstar
        wgt=readmask
        if Nobj>6:
            try:
                fmodel=stamp_name.replace("_star.fits","_modelxxxx.fits")
                imodel=fits.open(fmodel)[0].data
            except:
                print("generate wmodel",Nobj)
                print(FILTER,detector[DETECTOR],OPD)
                hdu=fits.open(stamp_name);print("hdu")
                imodel=np.zeros((Nobj,(Ng0+Nim)*osam,(Ng0+Nim)*osam),dtype='float')
                print("nrc")
                nrc = webbpsf.NIRCam();
                nrc.filter =FILTER
                nrc.detector=detector[DETECTOR]
                nrc.load_wss_opd_by_date(OPD)
                for ic in range(Nobj):
                    print(slicename+" mixpsf imodel",ic+1,spos[ic][0],spos[ic][1])
                    nrc.detector_position=(spos[ic][0],spos[ic][1])
                    psf = nrc.calc_psf(fov_pixels=Ng0+Nim,oversample=osam)
                    imodel[ic]=(psf[2].data/np.sum((psf[2].data))).T
                modelname=stamp_name.replace("_star.fits","_model.fits")
                print(modelname)
                _write_fits(modelname,imodel)
            modelname=stamp_name.replace("_star.fits","_model.fits")
            imodel=fits.open(modelname)[0].data
            PCs,pos,coeff,getNobj,slecstar,webcoeff,slecwgt,slecimodel=_web_psf_fit(
                instar,
                spos,
                NPCs,
                gain,
                method=1,
                wgt=wgt,
                imodel=imodel,
                osam=osam,
                SNRs=10.,
                big=big)
            PCname=stamp_name.replace("_star.fits","_PSF.fits")
            print(PCname,PCs.shape,pos.shape,coeff.shape,webcoeff.shape)
            _write_mult_fits(PCname,PCs,pos,coeff,webcoeff)
            hdu=fits.open(PCname)
            hdu[0].header=fits.open(stamp_name)[0].header
            hdu[0].header['OSAM']=osam
            hdu.writeto(PCname,overwrite=True)

            '''
            cent=np.zeros(2,dtype='float')
            chi=[]
            shape=[]
            R2=[]
            Ng0=instar.shape[1]
            sNg=60
            dpix=int((Ng0-sNg)/2);print("dpix:",dpix)
            rpsfs=_coeff2psf(pos,coeff,PCs,pos)
            residu=np.zeros((getNobj,sNg,sNg),dtype='float')
            for ic in range(getNobj):
                cutstar=slecstar[ic,dpix:dpix+sNg,dpix:dpix+sNg];print("slecstar")
                error_map=_wgtmap(cutstar);fsum=np.sum(cutstar);print("_wgtmap")
                cutstar/=np.sum(cutstar);error_map/=fsum;print("normalizing")
                cent[0],cent[1],sigma=_centriod_psf(cutstar);print('cent[0],cent[1]',cent[0],cent[1])
                oe1,oe2,or2=_size(cutstar,cent,3.5);print("shape:",oe1,oe2,or2)
                dx=cent[0]-(sNg/2.+0.5)+1;dy=cent[1]-(sNg/2.+0.5)+1;print('dx,dy=',dx,dy)
                tmp_mod=_interp_cubic(rpsfs[ic],sNg,dy,dx,osam=osam);print("_interp_cubic")
                tmodel=_interp_cubic(slecimodel[ic],sNg,dy,dx,osam=osam);
                tmp_mod=tmp_mod+tmodel
                cent[0],cent[1],sigma=_centriod_psf(tmp_mod);
                re1,re2,rr2=_size(tmp_mod,cent,3.5);
                chi.append(np.mean((cutstar-tmp_mod)**2/error_map**2))
                shape.append([oe1-re1,oe2-re2,(or2-rr2)/or2])
                R2.append(or2)
                print(oe1-re1,oe2-re2,(or2-rr2)/or2)
                cutstar/=np.sum(cutstar);tmp_mod/=np.sum(tmp_mod)
                residu[ic]=cutstar-tmp_mod
            shape=np.array(shape);chi2=np.array(chi);R2=np.array(R2)
            print('de1',np.mean(shape[:,0]),np.std(shape[:,0]))
            print('de2',np.mean(shape[:,1]),np.std(shape[:,1]))
            print('dr2',np.mean(shape[:,2]),np.std(shape[:,2]))
            print(np.mean(chi2))
            print("return")
            '''
            return(PCname)
        else:
            pass


def _get_psf(PSFhdu,ipos,psf_size,degrees=3):
    """
    Obtain the PSF model at target positions

    Parameters:
    -------
    PCs: numpy 3d array
        the principal components(PCs), which is the return from get_PC
    spos: numpy 2d array
        the selected star image positions, which is the return from get_PC
    coeff: numpy 2d array
        the coefficients and the corresponding errors of PCs, which is the return from get_PC
    ipos: numpy 2d array
        the target positions of the PSF model required
    psf_size: int
        the pixel size:psf_size*psf_size of psf model 
    degrees: int
        order of polynomial used to fit the PSF field variations

    Returns:
    ----------
    rPSF: numpy 3d array
        the constructed PSF model array, which are arrange as rPSF[index][x direction][y direction]
    """
    PCs,spos,coeff=PSFhdu[0].data,PSFhdu[1].data,PSFhdu[2].data
    rPSF=_coeff2psf(spos,coeff,PCs,ipos,degrees)
    return(rPSF)


def _psf_recon(spos,coeffs,PCs,webcoeff,slecimodel):
    """
    Obtain the PSF model at observed PSF star positions

    Parameters:
    -------
    spos: numpy 2d array
        PSF star position
    coeffs: numpy 2d array
        PC coeffs of stars
    PCs: numpy 3d array
        Principal components

    Returns:
    ----------
    rPSF: numpy 3d array
        the constructed PSF model array, which are arrange as rPSF[index][x direction][y direction]
    """
    iNstar=spos.shape[0];Nstar=spos.shape[0]
    Ng=PCs.shape[1]
    Nb=PCs.shape[0]
    icoeff=np.zeros((iNstar,Nb),dtype='float')
    for ic in range(Nb):
        icoeff[:,ic]=coeffs[:,ic,0];#print(icoeff[:,ic])
    rPSF=np.zeros((iNstar,Ng,Ng),dtype='float')
    for ic in range(iNstar):
        for ipc in range(Nb):
            rPSF[ic]+=(PCs[ipc]*icoeff[ic][ipc]+webcoeff[ic]*slecimodel[ic])
            #rPSF[ic]+=PCs[ipc]*coeffs[ic][ipc][0]
    return(rPSF)



def _coeff2psf(spos,coeffs,PCs,gpos,degrees):
    """
    interpolate the PSF models at target positions, which is provide for get_psf

    Parameters:
    -------
    spos: numpy 2d array
        the selected star image positions, which is the return from get_PC
    coeff: numpy 2d array
        the coefficients and the corresponding errors of PCs, which is the return from get_PC
    PCs: numpy 3d array
        the principal components(PCs), which is the return from get_PC
    gpos: numpy 2d array
        the target positions of the PSF model required
    degrees: int
        order of polynomial used to fit the PSF field variations

    Returns:
    ----------
    rPSF: numpy 3d array
        the constructed PSF model array, which are arrange as rPSF[index][x direction][y direction]
    """
    iNstar=gpos.shape[0];Nstar=spos.shape[0]
    Ng=PCs.shape[1]
    Nb=PCs.shape[0]
    icoeff=np.zeros((iNstar,Nb),dtype='float')
    for ic in range(Nb):
        icoeff[:,ic]=_polyfit2d(spos,coeffs[:,ic,0],coeffs[:,ic,1],gpos,degrees);#print(icoeff[:,ic])
    rPSF=np.zeros((iNstar,Ng,Ng),dtype='float')
    for ic in range(iNstar):
        for ipc in range(Nb):
            rPSF[ic]+=PCs[ipc]*icoeff[ic][ipc]
    return(rPSF)


def _polyfit2d(pos,val,error,tar_pos,polyorder):
    """
    constructing bivariate polynomials to fit the PSF field variations
    
    Parameters:
    -------
    pos: 2d numpy array
        star position
    val: 2d numpy array
        PC coefficients at star position
    error: 2d numpy array
        errors of val
    tar_pos: 2d numpy array
        target position
    polyorder: int
        order of polynomial used to fit the PSF field variations

    Returns:
    ----------
    ival: 2d numpy array
        PC coefficients at target position
    """
    Nstar=pos.shape[0];iNstar=tar_pos.shape[0]
    l0=0
    for i in range(polyorder+1):
        l0+=i+1
    #print("l0=",l0)
    xy=np.zeros((l0,Nstar),dtype='float')
    a=np.zeros((l0,l0),dtype='float')
    b=np.zeros(l0,dtype='float')
    ival=np.zeros(iNstar,dtype='float')
    ixy=np.zeros((l0,iNstar),dtype='float')
    for ic in range(Nstar):
        x=pos[ic][0];y=pos[ic][1]
        l=0
        for i in range(polyorder+1):
            for j in range(polyorder+1):
                if (i+j)<=polyorder:
                    xy[l][ic]=x**i*y**j 
                    l+=1
    for l in range(l0):
        for f in range(l0):
            for ic in range(Nstar):
                a[l][f]+=xy[l][ic]*xy[f][ic]*(1./error[ic]**2)
        for ic in range(Nstar):
            b[l]+=xy[l][ic]*val[ic]*(1./error[ic]**2)
    poly=np.linalg.solve(a,b);#print(poly)
    for ic in range(iNstar):
        x=tar_pos[ic][0];y=tar_pos[ic][1]
        l=0
        for i in range(polyorder+1):
            for j in range(polyorder+1):
                if (i+j)<=polyorder:
                    ixy[l][ic]=x**i*y**j 
                    l+=1
    for ic in range(iNstar):
        ival[ic]=np.sum(poly*ixy[:,ic])
    return(ival)


'''def coeff_map(file_PC,polyorder=1):
    nimgx,nimgy=
    pass
'''
def _PCA_star(instars,spos,NPCs,gain,psf_size,SNR_thresh):
    """
    constructing PSF model by using PCA

    Parameters:
    -------
    instars: 3d numpy array
        input PSF star images used to construct principal components
    spos: 2d numpy array
        position of PSF stars
    NPCs: int
        required number of principal components
    gain: float
        CCD gain
    psf_size: int 
        the pixel size:psf_size*psf_size of psf model 
    SNR_thresh: float
        signal to noise ratio threshold for the star images,which will used in PCs calculation


    Returns:
    ----------
    PCs: 3d numpy array
        the principal components(PCs)
    slecpos: 2d numpy array
        the selected star image positions aranged as spos[ic][0]->x,spos[ic][1]->y, 
        this will be used for the PSF field variation modelings
    coeff: 2d numpy array
        coefficients of PCs
    slecstar.shape[0]: int
        number of the selected star image
    slecstar: 3d numpy array
        slected star images
    """
    #sub background and algin center
    npc=NPCs
    Ng0=instars.shape[1]
    cent=np.zeros(2,dtype='float')
    pcastar=np.zeros((0,psf_size*psf_size),dtype='float')
    slecstar=np.zeros((0,Ng0,Ng0),dtype='float')
    tmp1=np.zeros((1,psf_size*psf_size),dtype='float')
    tmp2=np.zeros((1,Ng0,Ng0),dtype='float')
    ctr=np.zeros((0,psf_size,psf_size),dtype='float')
    tmp3=np.zeros((1,psf_size,psf_size),dtype='float')
    slecpos=[]
    ic=0;
    for star in instars:
        mean,sigma=_gaus_estimate(star*gain)
        star-=mean
        SNRs=_Singnal2Noise(star,gain=gain)
        if SNRs>SNR_thresh:
            # center algin
            star/=np.sum(star)
            cent[0],cent[1],sigma=_centriod_psf(star);#print('center',cent[0],cent[1])
            dx=cent[0]-(psf_size/2.+0.5);dy=cent[1]-(psf_size/2.+0.5);
            #print('dx,dy',dx,dy)
            #tmp_star=interp_cubic(star,psf_size,dy,dx,osam=1);#plt.imshow(tmp_star**0.1);plt.show()
            tmp_star=_Interp_bicubic(star,cent,psf_size,psf_size)
            tmp1[0]=tmp_star.ravel()
            pcastar=np.concatenate([pcastar,tmp1])
            tmp2[0]=star
            slecstar=np.concatenate([slecstar,tmp2])
            slecpos.append([spos[ic][1],spos[ic][2]])
            tmp3[0]=tmp_star
            ctr=np.concatenate([ctr,tmp3])
        ic+=1
    slecpos=np.array(slecpos)
    Nobj=slecstar.shape[0]
    coeff=np.ones((Nobj,npc,2),dtype='float')
    pcastar=pcastar.T
    PCs,sigma,v=np.linalg.svd(pcastar)
    #print(PCs.shape,sigma.shape,v.shape)
    for i in range(npc):
        coeff[:,i,0]=sigma[i]*v.T[:,i]
    PCs=(PCs[:,0:npc].T).reshape(npc,psf_size,psf_size)
    #write_fits("fits/star_ctr.fits",ctr)
    #print(PCs.shape,coeff.shape)
    return(PCs,slecpos,coeff,slecstar.shape[0],slecstar)

    
    



def _gaus_estimate(image):
    """
    estimate the background intensity and noise for the input star image
    Parameters:
    -------
    image: 2d numpy array
        image to estimate the noise

    Returns:
    ----------
    mean: float
        background intensity
    error: float 
        background noise
    """
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


def _Singnal2Noise(image,gain=1):
    """
    estimate the sognal to nosie ration of an image

    Parameters:
    -------
    image: 2d numpy array
        image
    gain: float
        CCD gain
    Returns:
    ----------
    snr: float
        SNR of image
    """
    Ng1=image.shape[0]
    Ng2=image.shape[1]
    weight=np.zeros(image.shape,dtype='float')
    mean,sigma=_gaus_estimate(image)
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



def _web_psf_fit(stars,spos,npc,gain,method,wgt=None,imodel=None,osam=1,SNRs=100.,big=0):
    '''
    in star image stamps: stars[Nstar][Ng][Ng]
    in star image positions: spos[Nstar][0]:x;spos[Nstar][1]:y
    number of PCs
    gain of CCD: gain
    method of reconstructions
    '''
    print("osam=%d"%(osam))
    Nstar=stars.shape[0];
    Ng0=stars.shape[1];
    dpix=4
    Nh=Ng0/2
    Ng=Ng0-dpix;
    Ngm=int((imodel.shape[1])/osam);print("Ngm",Ngm)
    print(Nstar,Ng0,npc,gain,method)
    fpoint=ctypes.POINTER(ctypes.c_float)
    ipoin1=ctypes.POINTER(ctypes.c_int)
    cint=ctypes.c_int
    cfloat=ctypes.c_float;
    _libpsffit.web_psf_fit.argtypes=[fpoint,fpoint,cint,cint,
    cint,cint,fpoint,fpoint,cfloat,cint,ipoin1,cint,cfloat,fpoint,fpoint,fpoint,cint]
    star=((Nstar)*Ng0*Ng0*ctypes.c_float)()
    mask=((Nstar)*Ng0*Ng0*ctypes.c_float)()
    PCs=(npc*Ng*osam*Ng*osam*ctypes.c_float)()
    rePCs=np.zeros((npc,Ng*osam,Ng*osam),dtype='float')
    starpos=(Nstar*2*ctypes.c_float)()
    webmodel=((Nstar)*Ngm*osam*Ngm*osam*ctypes.c_float)()
    wcoeff=(Nstar*2*ctypes.c_float)()
    
    Coeff=(Nstar*npc*2*ctypes.c_float)()
    
    Nobj=(1*ctypes.c_int)()
    print
    for ic in range(Nstar):
        k=ic*2
        starpos[k]=spos[ic][0]
        starpos[k+1]=spos[ic][1]
        for i in range(Ng0):
            for j in range(Ng0):
                k=ic*Ng0*Ng0+i*Ng0+j
                star[k]=stars[ic][i][j]


    if wgt is None:
        mask.all=1
    else:
        for ic in range(Nstar):
            for i in range(Ng0):
                for j in range(Ng0):
                    k=ic*Ng0*Ng0+i*Ng0+j
                    mask[k]=wgt[ic][i][j]

    if imodel is None:
        webmodel.all=0.
    else:
        for ic in range(Nstar):
            for i in range(Ngm*osam):
                for j in range(Ngm*osam):
                    k=ic*Ngm*osam*Ngm*osam+i*Ngm*osam+j
                    webmodel[k]=imodel[ic][i][j]
    print("start delivery")
    _libpsffit.web_psf_fit(star,starpos,npc,Nstar,Ng0,Ng,PCs,
        Coeff,gain,method,Nobj,osam,SNRs,webmodel,mask,wcoeff,big)

    #the side need to be changed
    reNobj=Nobj[0];print("reNobj:",reNobj,npc,Ng0,osam,Ngm)
    restarpos=np.zeros((reNobj,2),dtype="float")
    reCoeff=np.zeros((reNobj,npc,2),dtype="float")
    slecstar=np.zeros((reNobj,Ng0,Ng0),dtype="float")
    slecmask=np.zeros((reNobj,Ng0,Ng0),dtype="float")
    slecmodel=np.zeros((reNobj,Ngm*osam,Ngm*osam),dtype="float")
    webcoeff=np.zeros((reNobj,2),dtype="float")
    print("1")
    for ic in range(npc):
        for i in range(Ng*osam):
            for j in range(Ng*osam):
                k=ic*Ng*osam*Ng*osam+i*Ng*osam+j
                rePCs[ic][i][j]=PCs[k]
    print("2")
    for ic in range(reNobj):
        for i in range(npc):
            k=ic*npc*2+i*2
            reCoeff[ic][i][0]=Coeff[k]
            reCoeff[ic][i][1]=Coeff[k+1]
            #print(reCoeff[ic][i][0],)
        #print("\n")
        for i in range(2):
            k=ic*2
            restarpos[ic][0]=starpos[k];
            restarpos[ic][1]=starpos[k+1]
            webcoeff[ic][0]=wcoeff[k];
            webcoeff[ic][1]=wcoeff[k+1];
        for i in range(Ng0):
            for j in range(Ng0):
                slecstar[ic][i][j]=star[ic*Ng0*Ng0+i*Ng0+j]
                slecmask[ic][i][j]=mask[ic*Ng0*Ng0+i*Ng0+j]
        for i in range(Ngm*osam):
            for j in range(Ngm*osam):
                k=ic*Ngm*osam*Ngm*osam+i*Ngm*osam+j
                slecmodel[ic][i][j]=webmodel[k]
    
    print("reCoeff")
    return(rePCs,restarpos,reCoeff,reNobj,slecstar,webcoeff,slecmask,slecmodel)

def _sex_star(cat_path,cat_name,fits_path,fits_name,wdirs,stamp_size,magleft=None,magright=None,gain=None):
    """
    extract the star candidate from sextractor generated catalogue
    Parameters:
    -------
    cat_path: string
        path of sextractor generated catalogue
    cat_name: string
        name of sextractor generated catalogue
    fits_path: string 
        path of the fitsfile
    fits_name: string
        fitsfile name
    wdirs: string
        path to write the PSF star catalogue and images stamps
    stamp_size: string
        the pxiel size of the cutout star stamp images:stamp_size*stamp_size
    magleft: string
        a threshold to select the star images,which are defined log10(flux)
    magright: string
        a threshold to select the star images,which are defined log10(flux)
    gain: string

    Returns:
    ----------
    stamp_name：str
        the file name of star stamps
    cat_nam: str
        the catalogue name of star images
    """
    cint=ctypes.c_int
    cppchar = ctypes.POINTER(ctypes.POINTER(ctypes.c_char))
    _star_extrac.staridf.argtypes=[cint,cppchar]
    _star_extrac.staridf.restype = cint
    str_num=7
    if magleft is None:
        str_num=7
    else:
        if gain is None:
            str_num=9
        else:
            str_num=10
    
    
    strings = ["cut",cat_path,cat_name,fits_path,fits_name,wdirs,stamp_size,magleft,magright,gain]
    print(strings)
    c_strings = (str_num*ctypes.c_char_p)(*map(lambda s: s.encode(), strings))
    str_input = ctypes.cast(c_strings, cppchar)
    

    result=_star_extrac.staridf(str_num,str_input)






def _interp_cubic(image,target_Npix,dx,dy,osam):
    """
    shfit the image using cubic spline interpolation

    Parameters:
    -------
    image: 2d numpy array
        image to be shifted
    target_Npix: int 
        interpolated image pixel size
    dx: float
        shift distance in the first indx
    dy: float
        shift distance in the second indx
    osam: int
        ovsersampling factor of image

    Returns:
    ----------
    outI: 2d numpy array
        shifted image
    """
    #print(image.shape,target_Npix,dx,dy,osam)
    newf1=np.zeros(image.shape,dtype='float')
    Ng0=image.shape[0];Nge=target_Npix*osam
    dpix=((Ng0-Nge)/2);edx=dx*osam;edy=dy*osam
    outI=np.zeros((target_Npix,target_Npix),dtype='float')
    #print(dpix,edx,edy)
    for ic in range(Ng0):
        f1 = interpolate.interp1d(range(Ng0), image[ic] ,kind='cubic')
        newx=np.arange(0,Nge,1)+dpix-edx
        newf1[ic][0:Nge]=f1(newx)
    for ic in range(Nge):
        f2 = interpolate.interp1d(range(Ng0), newf1[:,ic], kind='cubic')
        newy=np.arange(0,Nge,1)+dpix-edy
        newf1[0:Nge,ic]=f2(newy)
    
    for i in range(target_Npix):
        for j in range(target_Npix):
            for it in range(osam):
                ra=int(i*osam+it);
                for jt in range(osam):
                    rb=int(j*osam+jt)
                    outI[i][j]+=newf1[ra][rb]
    return(outI)

def _Interp_bicubic(image,center,cutx,cuty):
    '''
    using the pixel center of image and cubic spline to interpolate the image (C code based)
    to the center of the pixel frame

    Parameters:
    -------
    image: 2d numpy array
        image to be shifted
    center: 1d numpy array
        center position in pixel frame
    cutx: int
        interpolated image pixel size in the 1st indx
    cuty: int
        interpolated image pixel size in the 2nd indx

    Returns:
    ----------
    out: 2d numpy array
        shifted image
    '''
    Ng1=image.shape[0]
    Ng2=image.shape[1]
    interpolated=(cutx*cutx*ctypes.c_float)()
    cent=(2*ctypes.c_double)(*list(np.ravel(center)))
    in_image=(Ng1*Ng2*ctypes.c_float)(*list(np.ravel(image)))
    _libpsffit.Interp_bicubiclib.argtypes=[fpoint,ctypes.c_int,ctypes.c_int,\
    ctypes.c_int,ctypes.c_int,dpoint,fpoint]
    _libpsffit.Interp_bicubiclib(in_image,Ng1,Ng2,cutx,cuty,cent,interpolated)
    out=np.array(list(interpolated)).reshape(cutx,cuty)[0:cutx]
    return(out)


def _centriod_psf(psf):
    '''
    intput the psf arry with Ng,Ng,
    and return the estimated pixel center of the profile :cx,cy
    based on fast shape estimate algorithm 

    Parameters:
    -------
    psf: 2d numpy array
        psf image

    Returns:
    ----------
    cent[0]: float
        center position of 1st array indx
    cent[1]: float
        center position of 2nd array indx
    cent[2]: float
        gaussian fitted sigma of psf
    '''
    Ng1=psf.shape[0]
    Ng2=psf.shape[1];
    in_image=(Ng1*Ng2*ctypes.c_float)(*list(np.ravel(psf)))
    cent=(3*ctypes.c_double)()
    _libpsffit.centriod_psflib.argtypes=[fpoint,ctypes.c_int,ctypes.c_int,dpoint]
    _libpsffit.centriod_psflib(in_image,Ng1,Ng2,cent)
    return(cent[0],cent[1],cent[2])



def _size(image,center,sigma):
    '''
    estimated the shape parameters of image by using the second brightness moment

    Parameters:
    -------
    image:2d numpy array
        image to estimated
    center: 1d numpy array
        center position of image, which is estimated from centriod_psf()
    sigma: float
        width of the gaussian weight function

    Returns:
    ----------
    e1: float
        the 1st ellipticity component
    e2: float
        the 2nd ellipticity component
    R2: float
        the size of image 
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
        print("R2<=0\n");return(e1,e2,r2)
    else:
        return(e1,e2,R2)



def _write_fits(fitsname,data,head1=None,head2=None):
    """
    write data in fits file

    Parameters:
    -------
    fitsname: str 
        name of fits file
    data: numpy array
        data of fits file
    head1: fits header
        header of fisrt fits hdu list
    head2: header
        header of fisrt fits hdu list

    Returns:
    ----------
    """
    wdata=data;
    if head1==None :
        #print("head1 none",fitsname)
        hdu = fits.PrimaryHDU(wdata)
        hdul = fits.HDUList([hdu])
        hdul.writeto(fitsname,overwrite=True)
    if head1!=None and head1==None:
        #print("head2 none",fitsname)
        hdu = fits.PrimaryHDU(wdata,header=head1)
        hdul = fits.HDUList([hdu])
        hdul.writeto(fitsname,overwrite=True)
    if head1!=None and head2!=None:
        #print("head1 not none",fitsname)
        hdu = fits.PrimaryHDU(header=head1)
        hdu1= fits.ImageHDU(wdata,header=head2)
        hdul = fits.HDUList([hdu,hdu1])
        hdul.writeto(fitsname,overwrite=True)


def _write_mult_fits(fitsname,data1=[],data2=[],data3=[],data4=[],data5=[]):
    """
    write multiple array data in one fits file

    Parameters:
    -------
    fitsname: str 
        name of fits file
    data1: numpy array
        data of fits file
    data2: numpy array
        data of fits file
    data3: numpy array
        data of fits file

    Returns:
    ----------

    """
    if type(data2) is list :
        hdu = fits.PrimaryHDU(data1)
        hdul = fits.HDUList([hdu])
        hdul.writeto(fitsname,overwrite=True) 
    if type(data2) is np.ndarray and type(data3) is list:
        hdu = fits.PrimaryHDU(data1)
        hdu1=fits.ImageHDU(data2)
        hdul = fits.HDUList([hdu,hdu1])
        hdul.writeto(fitsname,overwrite=True) 
    if type(data3) is np.ndarray and type(data4) is list:
        hdu = fits.PrimaryHDU(data1)
        hdu1=fits.ImageHDU(data2)
        hdu2=fits.ImageHDU(data3)
        hdul = fits.HDUList([hdu,hdu1,hdu2])
        hdul.writeto(fitsname,overwrite=True)
    if type(data4) is np.ndarray and type(data5) is list:
        hdu = fits.PrimaryHDU(data1)
        hdu1=fits.ImageHDU(data2)
        hdu2=fits.ImageHDU(data3)
        hdu3=fits.ImageHDU(data4)
        hdul = fits.HDUList([hdu,hdu1,hdu2,hdu3])
        hdul.writeto(fitsname,overwrite=True)
    if type(data5) is np.ndarray :
        hdu = fits.PrimaryHDU(data1)
        hdu1=fits.ImageHDU(data2)
        hdu2=fits.ImageHDU(data3)
        hdu3=fits.ImageHDU(data4)
        hdu4=fits.ImageHDU(data5)
        hdul = fits.HDUList([hdu,hdu1,hdu2,hdu3,hdu4])
        hdul.writeto(fitsname,overwrite=True)


def _wgtmap(image,gain=1):
    """
    estimate the error image of image

    Parameters:
    -------
    image: 2d numpy array
        image to be estimated
    gain: float
        CCD gain

    Returns:
    ----------
    weight: 2d numpy array
        error image of input image
    """
    Ng1=image.shape[0]
    Ng2=image.shape[1]
    weight=np.zeros(image.shape,dtype='float')
    mean,sigma=_gaus_estimate(image)
    #print("mean=%f,sigma=%f"%(mean,sigma))
    for i in range(Ng1):
        for j in range(Ng2):
            detx=sigma*sigma*gain*gain
            dety=gain*np.fabs(image[i][j])
            if dety<detx : error=sigma
            else : error=np.sqrt(dety+detx)/gain
            #weight[i][j]=1./(error**2)
            weight[i][j]=error
    return(weight)


def _S2N(image,gain=1):
    Ng1=image.shape[0]
    Ng2=image.shape[1]
    weight=np.zeros(image.shape,dtype='float')
    mean,sigma=_gaus_estimate(image)
    image-=mean
    #mean=sigma=0
    #print("method=%f,sigma=%f"%(mean,sigma))
    for i in range(Ng1):
        for j in range(Ng2):
            detx=sigma*sigma*gain*gain
            dety=gain*np.fabs(image[i][j])
            if dety<detx : error=sigma
            else : error=np.sqrt(dety+detx)/gain
            weight[i][j]=error
    snr=np.sum(image)/np.sqrt(np.sum(weight**2))
    return(snr)

def _get_webPSF(PCfile,
              psfpos,
              Target_Ng):
    '''
    Using the calculated PCs to interpolate the PSF modeling 

    parameters
    ----------
    PCfile: str
        the constructed PC file
    psfpos: np.ndarray
        the required psf pos which you want, pos can be a 2 dimentional array,such as psfpos[xpos][ypos]
    Target_Ng: int
          the required pixel size of psf stamp

    return
    ------
    np.ndarray:
        the interpolated PSF modeling image matrix using polynomials
    '''
    print("in get_webPSF")
    rNstar=6
    tmpname=PCfile.split("/");l0=len(tmpname);l1=len(tmpname[l0-1]);
    slicename=tmpname[l0-1][0:l1-len(".fits")]
    try:
        hdu=fits.open(PCfile);
    except:
        print("file: "+PCfile+" open failed, please check the file name!")
        sys.exit()
    current_dir=os.getcwd()
    Nobj=psfpos.shape[0];psfcat=np.zeros((Nobj,2));print("Nobj:",Nobj)
    for ic in range(Nobj):
        psfcat[ic][0]=psfpos[ic][0];psfcat[ic][1]=psfpos[ic][1];
    if os.path.exists(PCfile)==False :
        print("file does not exists");
        pass
    else:
        # use mixture psf model
        print("use mixture psf model")
        OPDs=hdu[0].header['DATE-BEG']
        Filter=hdu[0].header['FILTER']
        Detect=hdu[0].header['DETECTOR']
        Osam=hdu[0].header['OSAM'];print(Filter,Detect,OPDs,Osam)
        PCs=hdu[0].data;
        npc=PCs.shape[0]
        spos=hdu[1].data
        coeff=hdu[2].data
        webcoeff=hdu[3].data
        Nstar=spos.shape[0];
        print(PCs.shape,npc,spos.shape,coeff.shape,webcoeff.shape,Target_Ng*Osam)
        rPSF=np.zeros((Nobj,int(Target_Ng*Osam),int(Target_Ng*Osam)),dtype='float')
        print("make rPSF")
        rPSF=np.zeros((Nobj,int(Target_Ng*Osam),int(Target_Ng*Osam)),dtype='float');
        print("rPSF")
        Tmodel=np.zeros((Nobj,Target_Ng,Target_Ng),dtype='float');print("Tmodel")
        print(Nstar,rNstar)
        if Nstar<rNstar:
            print("too few star, use webbpsf")
            pass
        else:
            print("hybird psf model")
            #coeff=hdu[2].data;print(psfpos.shape)
            polyorder=2
            if Nstar>10:polyorder=3
            #if Nstar>15:polyorder=4
            #if Nstar>21:polyorder=5
            print("_coeff2psf")
            print(spos.shape,coeff.shape,PCs.shape,psfpos.shape,polyorder)
            rPSF=_coeff2psf(spos,coeff,PCs,psfpos,degrees=polyorder);print("rPSF",rPSF.shape)
            eNg=rPSF.shape[1];print("eNg=",eNg)
            sNg=eNg;dNg=int((eNg-sNg)/2)
            cut_rPSF=np.zeros(rPSF.shape,dtype='float')
            for ic in range(Nobj):
                for i in range(eNg):
                    ra=i-(eNg/2-0.5)
                    for j in range(eNg):
                        rb=j-(eNg/2-0.5)
                        r2=ra**2+rb**2
                        r=r2**0.5
                        if r2<(sNg/2)**2:
                            cut_rPSF[ic][i][j]=1
                            cut_rPSF[ic][i][j]=1./(1+np.exp((3.*(r-9))))
                cut_rPSF[ic]*=rPSF[ic]
            icoeff=_polyfit2d(spos,webcoeff[:,0],webcoeff[:,1],psfpos,polyorder=3)
            print(icoeff)
            if Target_Ng<=eNg/Osam:
                print("Target_Ng:%d\t%d"%(Target_Ng,eNg/Osam));
                cent=np.zeros(2,dtype='float');
                print(Filter,Detect,OPDs,Osam)
                webp=webbpsf.NIRCam()
                webp.filter=Filter
                webp.detector=detector[Detect]
                webp.load_wss_opd_by_date(OPDs)
                dpix=int((eNg-Target_Ng*Osam)/2);print("dpix:",dpix)
                for ic in range(Nobj):
                    print("ic:",ic,Nobj,psfcat[ic][0],psfcat[ic][1])
                    webp.detector_position=(psfcat[ic][0],psfcat[ic][1]);print("position")
                    psf = webp.calc_psf(fov_pixels=int(eNg/Osam),oversample=Osam);print("calc_psf")
                    psf=psf[2].data.T;psf/=np.sum(psf);print("normalizing")
                    #psf*=0
                    '''plt.imshow((np.fabs(cut_rPSF[ic,:,:]))**0.1);plt.show()
                    plt.imshow((np.fabs(psf[:,:]*icoeff[ic]))**0.1);plt.show()'''
                    tmppsf=psf[:,:]*icoeff[ic]+cut_rPSF[ic,:,:];print("addition")
                    '''
                    for i in range(eNg):
                        for j in range(eNg):
                            if  tmppsf[i][j]<0 : #padding the nagetive value
                                tmppsf[i][j]=0
                    '''
                    indx=np.where(tmppsf<0.)
                    tmppsf[indx]=0
                    for i in range(eNg):
                        for j in range(eNg):
                            if  tmppsf[i][j]==0 and i-1>0 and i+1<eNg and j-1>0 and j+1<eNg: #padding the nagetive value
                                tmppsf[i][j]=(tmppsf[i-1][j-1]+tmppsf[i-1][j]+tmppsf[i-1][j+1]+
                                    tmppsf[i][j-1]+tmppsf[i][j+1]+
                                    tmppsf[i+1][j-1]+tmppsf[i+1][j]+tmppsf[i+1][j+1])/8.
                    tmppsf/=np.sum(tmppsf);print("padding")


                    Tmodel[ic]*=0
                    tmp=tmppsf[dpix:dpix+Target_Ng*Osam,dpix:dpix+Target_Ng*Osam]
                    for i in range(Target_Ng):
                        for j in range(Target_Ng):
                            for l in range(Osam):
                                ra=i*Osam+l;
                                for m in range(2):
                                    rb=j*Osam+m
                                    Tmodel[ic][i][j]+=tmp[ra][rb]
                    #print(Tmodel[ic][125][125])
            else:  #Target_Ng>=eNg/osam
                print(Target_Ng,eNg)
                dpix=int((Target_Ng*Osam-eNg)/2);print("case 2")
                webp=webbpsf.NIRCam()
                webp.filter=Filter
                webp.detector=detector[Detect]
                webp.load_wss_opd_by_date(OPDs)
                print("dpix:",dpix)
                for ic in range(Nobj):
                    webp.detector_position=(psfcat[ic][0],psfcat[ic][1]) 
                    psf = webp.calc_psf(fov_pixels=Target_Ng,oversample=Osam)
                    psf=psf[2].data.T;psf/=np.sum(psf)
                    psf/=np.sum(psf[dpix:dpix+eNg,dpix:dpix+eNg]) #calibrate the flux ratio
                    psf*=icoeff[ic]
                    psf[dpix:dpix+eNg,dpix:dpix+eNg]+=cut_rPSF[ic,:,:]
                    tmp=psf
                    '''
                    for i in range(eNg):
                        for j in range(eNg):
                            if  tmp[i+dpix][j+dpix]<0 : #padding the nagetive value
                                tmp[i+dpix][j+dpix]=0
                    '''
                    indx=np.where(tmp<0.)
                    if indx[0].shape[0]>0:
                        tmp[indx]=0
                    #plt.imshow(((tmp))**0.1);plt.show()
                    for i in range(Target_Ng*Osam):
                        for j in range(Target_Ng*Osam):
                            if  tmp[i][j]==0 and i-1>0 and i+1<Target_Ng*Osam and j-1>0 and j+1<Target_Ng*Osam: #padding the nagetive value
                                tmp[i][j]=(tmp[i-1][j-1]+tmp[i-1][j]+tmp[i-1][j+1]+tmp[i][j-1]+tmp[i][j+1]+tmp[i+1][j-1]+tmp[i+1][j]+tmp[i+1][j+1])/8.
                    tmp/=np.sum(tmp)
                    Tmodel[ic]*=0
                    for i in range(Target_Ng):
                        for j in range(Target_Ng):
                            for l in range(Osam):
                                ra=i*Osam+l;
                                for m in range(2):
                                    rb=j*Osam+m
                                    Tmodel[ic][i][j]+=tmp[ra][rb]
                    Tmodel[ic]/=np.sum(Tmodel[ic]) 
    return(Tmodel)

