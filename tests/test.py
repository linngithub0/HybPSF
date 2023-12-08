from HybPSF.HybPSF import _Image2PSFmodel,_get_webPSF
import numpy as np
import matplotlib.pyplot as plt
'''fname="/Users/linn/Documents/soft/HybPSF/tests/jw02736001001_02105_00009_nrcalong_cal.fits"
wdirs="/Users/linn/Documents/soft/HybPSF/tests/"
PCname=_Image2PSFmodel(fname,wdirs,NPCs=5,SNR_thresh=10.,osam=4)
print((PCname))'''
PCname="/Users/linn/Documents/soft/HybPSF/tests/catalogue/star_stamps/jw02736001001_02105_00009_nrcalong_cal_PSF.fits"
pos=np.array([[0,0]])
rPSFs=_get_webPSF(PCname,pos,26);
print(rPSFs.shape)
plt.imshow(rPSFs[0])
plt.show()
