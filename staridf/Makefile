#cfitsio library path, library name
FITSPATH=-L/home/lnie/code/jwst_psf/lib
FITSLIB=-lcfitsio
FITSINCLUDE=-I/home/lnie/code/jwst_psf/include
# general include path
INCLUDE=-I/usr/local/include/
STDLIB=-L/usr/local/lib/
 
objects = main.o fits_IO.o nrutil.o wavelet_Trous.o sort.o path_option.o SNRs.o fof_source.o read_funcpsfs.o ran1.o gasdev.o
  
cut:  $(objects)
	cc -O3 -o $@ $(objects) $(FITSPATH) $(FITSLIB) -lm
# general compilation rules 
.SUFFIXES:  .c   .o
.c.o:
	cc -Wno-nullability-completeness -c $< -O3 -o $@ $(FITSINCLUDE) $(INCLUDE)
.PHONY:clean
clean:
	-rm cut
	-rm ${objects}
