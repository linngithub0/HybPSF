#configurations for the star catalogue
#file_name: the fits file
file_name="./test/jw02736001001_02103_00001_nrcalong_cal_1overf.fits"

#cata_path:  output path used to write the catalogue files
cata_path="/home/lnie/code/jwst_psf/test/"

#file_name="/data/JWST_rev/NIRCam/simp_coadd/fiDrizzleRC/F090Wb2_fiDrizzleRC.fits"

#npc:  number of PCs, for current method, 4 are recommend
npc=4

#diagnose:  if diagnose==True, some diagnose file will be written in cata_path
diagnose=True

#osam:  the oversampleing factor, 1 or 2 could be selected
osam=2

#the following are the instrument configurations
Instr="NIRCAM"
Detec="NRCBLONG"
Filt="F444W"




