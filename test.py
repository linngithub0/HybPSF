from HybPSF.HybPSF import _Image2PSFmodel,_get_psf
import numpy as np

fname="/Users/linn/Documents/soft/HybPSF/tests/jw02736001001_02105_00009_nrcalong_cal.fits"
wdirs="/Users/linn/Documents/soft/HybPSF/tests/"
PSFhdu=_Image2PSFmodel(fname,wdirs,NPCs=5,10.,osam=4)
print(type(PSFhdu))
rPSFs=_get_psf(PSFhdu,np.array([[200,100]]),psf_size=25)
plt.imshow(rPSFs[0])
plt.show()