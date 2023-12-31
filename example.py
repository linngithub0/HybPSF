#for JWST
import numpy as np
import matplotlib.pyplot as plt
from psf_fit import get_webPSF,write_fits
#the position of the PSF should be arranged into a numpy array
cat=np.array([[100,300],[500,320]])
#the file path of the corresponding fits, which are used to provide the detector information
#the PC file path, which are located at the followings, you can also download them to your own device, and then change this path
PCfile="./test/jw02736001001_02103_00001_nrcalong_cal_1overf_PC.fits"
#use get_webPSF to get the returned PSF model, if the PC file of the fitsfile included in the PC path, which is caused by the too few number of stars, the returned PSF are generated by WebbPSF
rPSF=get_webPSF(PCfile=PCfile,
              psfpos=cat,
              Target_Ng=200,
              osam=1)
for stamp in rPSF:
    plt.imshow(np.log(stamp))
    plt.show()
