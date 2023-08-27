#cfitsio library path, library name
FITSPATH=-L/home/lnie/code/jwst_psf/lib
FITSLIB=-lcfitsio
FITSINCLUDE=-I/home/lnie/code/jwst_psf/include 
# general include path
INCLUDE=-I/usr/local/include/
STDLIB=-L/usr/local/lib/
 
objects=mt.o read_funcpsfs_new.o fits_IO.o nrutil.o creat_basisfs.o estimate_para.o fit_gaus.o gasdev.o gaus_estimate.o indexx.o iSPCA_entr.o itoa.o k_sigma.o lubksb.o ludcmp.o moffat.o moffatlets.o mov_mer_ba.o PCA.o method_branches.o iSPCA.o SPCA_entr.o C_tools_lib.o ran1.o svdcmp.o rebasis.o star_ellip.o pythag.o linear_operation.o EMPCA.o EMPCA_entr.o interp_bicubic.o SPCA.o npsf.o coeff_interp.o size.o sort.o creat_star.o krig.o  fof_source.o
  
libpsffit.so:  $(objects) 
	cc -fPIC -g -O3 -o $@ $(objects) $(FITSPATH) $(FITSLIB) $(STDLIB) -shared -lm
 
.SUFFIXES:  .c   .o
.c.o:
	cc -fPIC -g -Wno-nullability-completeness -c $< -O3 -o $@ $(FITSINCLUDE) $(INCLUDE)
.PHONY:clean
clean:
	-rm libpsffit.so
	-rm ${objects}
