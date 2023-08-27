rm *.so
rm *.dylib
rm -r lib
rm -r include
cd source
make clean
cd ..
a=`pwd`
cd fftw-3.3.10/
./configure --prefix=$a\/
make clean
make
make install
cd ../cfitsio-4.1.0/
./configure --prefix=$a\/
make clean
make
make install
cd ..
dir=`pwd`
cd source
echo "#cfitsio library path, library name" > Makefile
echo "FITSPATH=-L"${dir}"/lib" >> Makefile
echo "FITSLIB=-lcfitsio -lcurl" >> Makefile
echo "FITSINCLUDE=-I"${dir}"/include" >> Makefile
echo "# general include path" >> Makefile
echo "INCLUDE=-I/usr/local/include/" >> Makefile
echo "STDLIB=-L/usr/local/lib/" >> Makefile
echo " " >> Makefile
echo "objects=mt.o read_funcpsfs_new.o fits_IO.o nrutil.o creat_basisfs.o estimate_para.o fit_gaus.o\
 gasdev.o gaus_estimate.o indexx.o iSPCA_entr.o itoa.o k_sigma.o lubksb.o ludcmp.o moffat.o\
 moffatlets.o mov_mer_ba.o PCA.o method_branches.o iSPCA.o SPCA_entr.o C_tools_lib.o\
 ran1.o svdcmp.o rebasis.o star_ellip.o pythag.o linear_operation.o EMPCA.o EMPCA_entr.o\
 interp_bicubic.o SPCA.o npsf.o coeff_interp.o size.o sort.o creat_star.o krig.o " >> Makefile
echo "  " >> Makefile
echo "libpsffit.so:  \$(objects) " >> Makefile
echo "	cc -fPIC -g -O3 -o \$@ \$(objects) \$(FITSPATH) \$(FITSLIB) \$(STDLIB) -shared -lm" >> Makefile
echo " " >> Makefile
echo ".SUFFIXES:  .c   .o" >> Makefile
echo ".c.o:" >> Makefile
echo "	cc -fPIC -g -Wno-nullability-completeness -c \$< -O3 -o \$@ \$(FITSINCLUDE) \$(INCLUDE)" >> Makefile
echo ".PHONY:clean" >> Makefile
echo "clean:
	-rm libpsffit.so
	-rm \${objects}" >> Makefile
make clean
make
cp libpsffit.so ../
cd ..
cp lib/libcfitsio.9.dylib ./



cd source/staridf/
echo "#cfitsio library path, library name" > Makefile
echo "FITSPATH=-L"${dir}"/lib" >> Makefile
echo "FITSLIB=-lcfitsio -lcurl" >> Makefile
echo "FITSINCLUDE=-I"${dir}"/include" >> Makefile
echo "# general include path" >> Makefile
echo "INCLUDE=-I/usr/local/include/" >> Makefile
echo "STDLIB=-L/usr/local/lib/" >> Makefile
echo " " >> Makefile
echo "objects = main.o fits_IO.o nrutil.o wavelet_Trous.o sort.o \
path_option.o SNRs.o fof_source.o read_funcpsfs.o ran1.o gasdev.o" >> Makefile
echo "  " >> Makefile
echo "cut:  \$(objects)" >> Makefile
echo "	cc -O3 -o \$@ \$(objects) \$(FITSPATH) \$(FITSLIB) -lm" >> Makefile
echo "# general compilation rules " >> Makefile
echo ".SUFFIXES:  .c   .o" >> Makefile
echo ".c.o:" >> Makefile
echo "	cc -Wno-nullability-completeness -c \$< -O3 -o \$@ \$(FITSINCLUDE) \$(INCLUDE)" >> Makefile
echo ".PHONY:clean" >> Makefile
echo "clean:
	-rm cut
	-rm \${objects}" >> Makefile
make clean
make
cp cut ../../










