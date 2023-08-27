import numpy as np
from psf_fit import get_webPSF,write_fits,poly_coeff_map
'''cat=np.array([100,300])
Tmodel=get_webPSF(Cam="NIRCam",
    Filter="F090W",
    detector="nrcb4",
    expos=1,
    psfpos=cat,
    Target_Ng=31,
    osam=2,
    sign="single")
write_fits("test.fits",Tmodel)
'''
'''decoeff=poly_coeff_map(Cam="NIRCam",
    Filter="F090W",
    detector="nrcb4",
    expos=1,
    polyorder=1)
print(decoeff.shape)
print(decoeff)
'''