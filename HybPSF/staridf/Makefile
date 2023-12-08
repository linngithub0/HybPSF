#cfitsio library path, library name
FITSPATH=-L/Users/linn/Documents/soft/HybPSF/HybPSF/lib
FITSLIB=-lcfitsio -lcurl
FITSINCLUDE=-I/Users/linn/Documents/soft/HybPSF/HybPSF/include
# general include path
INCLUDE=-I/usr/local/include/
STDLIB=-L/usr/local/lib/
 
objects = staridf.o fits_IO.o nrutil.o wavelet_Trous.o sort.o path_option.o SNRs.o fof_source.o read_funcpsfs.o ran1.o gasdev.o
  
star_extrac.so:  $(objects)
	cc -fPIC -g -O3 -o $@ $(objects) $(FITSPATH) $(FITSLIB) -shared -lm
# general compilation rules 
.SUFFIXES:  .c   .o
.c.o:
	cc -fPIC -g -Wno-nullability-completeness -c $< -O3 -o $@ $(FITSINCLUDE) $(INCLUDE)
.PHONY:clean
clean:
	-rm star_extrac.so
	-rm ${objects}
