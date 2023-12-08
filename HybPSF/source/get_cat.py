import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os
import time
import webbpsf
nrca_short_filters = ['F070W', 'F090W', 'F115W', 'F140M', 'F150W2', 'F150W', 'F162M', 'F164N', 'F182M', 'F187N',
                          'F200W', 'F210M', 'F212N']
nrca_long_filters = ['F250M', 'F277W', 'F300M', 'F322W2', 'F323N', 'F335M', 'F356W', 'F360M', 'F405N', 'F410M',
                         'F430M', 'F444W', 'F460M', 'F466N', 'F470N', 'F480M']

nrca_short_detectors = ['NRCA1', 'NRCA2', 'NRCA3', 'NRCA4', 'NRCB1', 'NRCB2', 'NRCB3', 'NRCB4']
nrca_long_detectors = ['NRCA5', 'NRCB5']
F090W_fname="jw02736001001_02101_"
F150W_fname="jw02736001001_02103_"
F200W_fname="jw02736001001_02105_"
F277W_fname="jw02736001001_02101_"
F356W_fname="jw02736001001_02103_"
F444W_fname="jw02736001001_02105_"
filte_name=dict({"F090W/":"jw02736001001_02101_","F150W/":"jw02736001001_02103_",
                "F200W/":"jw02736001001_02105_","F277/":"jw02736001001_02101_",
                "F356/":"jw02736001001_02103_","F444W/":"jw02736001001_02105_"})
detector=dict({"nrca1":"nrca1","nrca2":"nrca2","nrca3":"nrca3","nrca4":"nrca4",
               "nrcb1":"nrcb1","nrcb2":"nrcb2","nrcb3":"nrcb3","nrcb4":"nrcb4",
               "nrcalong":"nrca5","nrcblong":"nrcb5"})
gains=dict({"nrca1":"2.08","nrca2":"2.02","nrca3":"2.17","nrca4":"2.02",
            "nrcb1":"2.01","nrcb2":"2.14","nrcb3":"1.94","nrcb4":"2.03",
            "nrcalong":"1.84","nrcblong":"1.80"})
dirs="data/JWST/"
Cam="NIRCam/"
filte="F090W/"
#detect=["nrcalong","nrcblong"]
detect=["nrca1","nrca2","nrca3","nrca4","nrcb1","nrcb2","nrcb3","nrcb4"]
append="_cal_1overf.fits"  #nrca4_cal_1overf.fits
wdirs="data/JWST/"
cat_header=["#   1 NUMBER                 Running object number\n",                                  
"#   2 X_IMAGE                Object position along x                                      [pixel]\n",
"#   3 Y_IMAGE                Object position along y                                      [pixel]\n",
"#   4 FWHM_IMAGE             FWHM assuming a gaussian core                                [pixel]\n",
"#   5 FLUX_RADIUS            Fraction-of-light radii                                      [pixel]\n",
"#   6 FLUX_AUTO              Flux within a Kron-like elliptical aperture                  [count]\n",
"#   7 FLUXERR_AUTO           RMS error for AUTO flux                                      [count]\n",
"#   8 FLAGS                  Extraction flags"                                          
"#   9 CLASS_STAR             S/G classifier output\n"]
try:
    cmd="mkdir /home/lnie/data/JWST_1overf/"
    os.system(cmd)
    cmd="mkdir /home/lnie/data/JWST_1overf/"+Cam
    os.system(cmd)
    cmd="mkdir /home/lnie/data/JWST_1overf/"+Cam+filte
    os.system(cmd)
    for det in detect:
        cmd="mkdir /home/lnie/data/JWST_1overf/"+Cam+filte+det
        os.system(cmd)
except:
    pass

#defination useful functions
def get_filename(file_dir,append):
    fname=[];
    for list in os.listdir(file_dir):
        if os.path.isfile(os.path.join(file_dir,list)):
            if list[(len(list)-len(append)):len(list)]==append:
                fname.append(list);

    return fname
def write_fits(fitsname,data,head1=None,head2=None):
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
def fft_conv(image,psf):
    '''
    convlove image with a psf,usually image size should be larger than psf
    '''
    size=image.shape
    core=np.zeros(size,dtype='float')
    nx=image.shape[0];ny=image.shape[1]
    px=psf.shape[0];py=psf.shape[1]
    ra=int((nx-px)/2);rb=int((ny-py)/2)
    core[ra-1:ra+px-1,rb-1:rb+py-1]=psf
    Fg=np.fft.fft2(image)
    Fp=np.fft.fft2(core)
    Fg*=Fp
    for i in range(Fp.real.shape[0]):
        for j in range(Fp.real.shape[1]):
            l=(i+j+2)%2
            if (l==1) :Fg.real[i][j]*=-1;Fg.imag[i][j]*=-1;
    tmp=np.fft.ifft2(Fg)
    cimage=tmp.real
    return(cimage)

#this part is to construct wgt map and get calaogue from sextractor
dirs="data/JWST/"
wdirs="data/JWST/"
append="_00007_nrca2_cal_1overf.fits"
for det in detect:
    detd=det+"/"
    fdirs=dirs            #contruct file directory
    print(fdirs)
    getname=get_filename(fdirs,append)                  #get the compelet file name with appendix "_cal.fits" in directory fdirs
    for gname in getname:
        tmpname=gname.split("/");l0=len(tmpname);l1=len(tmpname[l0-1]) #get only the filename of fits without directory
        slicename=tmpname[l0-1][0:l1-len(".fits")]
        fname=fdirs+gname
        print("fname:",fname) 
        
        hdu=fits.open(fname)                 #read fits
        data=hdu[1].data                     #get fits data
        wname=wdirs+getname[0]
        write_fits(wname,data,head1=hdu[1].header,head2=hdu[1].header)
        
        wgt=np.ones(data.shape,dtype='float')
        wgtname=wdirs+slicename+"_wgt.fits"
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                if data[i][j]<=0:
                    wgt[i][j]=0
        write_fits(wgtname,wgt)
        
        coldcat=wdirs+"cold.cat";print(coldcat)
        segname=wdirs+slicename+"_seg.fits"
        coldcmd="sex "+wname+" -CATALOG_NAME "+coldcat+" -DETECT_MINAREA 120\
        -DETECT_THRESH 3.2 -DEBLEND_NTHRESH 32 -DEBLEND_MINCONT 0.04 -BACK_SIZE 400\
        -BACK_FILTERSIZE 5 -BACKPHOTO_TYPE LOCAL -CHECKIMAGE_TYPE SEGMENTATION -CHECKIMAGE_NAME "+segname
        print("coldcmd",coldcmd)
        os.system(coldcmd)
        hotcat=wdirs+"hot.cat"
        hotcmd="sex "+wname+" -CATALOG_NAME "+hotcat+" -DETECT_MINAREA 18 -DETECT_THRESH 1.2\
        -DEBLEND_NTHRESH 32 -DEBLEND_MINCONT 0.065 -BACK_SIZE 100 -BACK_FILTERSIZE 3 -BACKPHOTO_TYPE LOCAL"
        os.system(hotcmd)
        try:
            print("start erosion")
            erosion_core=np.ones((40,40))
            hdu=fits.open(segname)
            segmap=hdu[0].data
            erosion=fft_conv(segmap,erosion_core)
            boxname=wdirs+slicename+"cat.pdf"
            cold_cat=np.loadtxt(coldcat);print("cold:",cold_cat.shape)
            hot_cat=np.loadtxt(hotcat);print("hot:",hot_cat.shape)
            Nobj=hot_cat.shape[0]
            hcopy=np.zeros(hot_cat.shape,dtype='float')
            k=0
            for it in range(Nobj):
                x=int((hot_cat[it][1]));y=int((hot_cat[it][2]))
                if erosion[x][y]<1.0e-5 :
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
            r50=cold_cat[:,4];fwhm=cold_cat[:,3];mag=-2.5*np.log10(cold_cat[:,5])
            fig=plt.figure()
            plt.plot(mag,r50,'.',c='blue',alpha=0.3)
            plt.plot(mag,fwhm,'.',c='red',alpha=0.5)
            plt.ylim(-2,30)
            plt.savefig(wdirs+slicename+'_cat.pdf')
            
        except:
            print("file open failed, please check this file %s"%(fname))
            pass
      