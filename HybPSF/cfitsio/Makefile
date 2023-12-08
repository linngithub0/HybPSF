#
# Makefile for cfitsio library:
#
# Oct-96 : original version by 
#
#       JDD/WDP
#       NASA GSFC
#       Oct 1996
#
# 25-Jan-01 : removed conditional drvrsmem.c compilation because this
#             is now handled within the source file itself.
# 09-Mar-98 : modified to conditionally compile drvrsmem.c. Also
# changes to target all (deleted clean), added DEFS, LIBS, added
# DEFS to .c.o, added SOURCES_SHMEM and MY_SHMEM, expanded getcol*
# and putcol* in SOURCES, modified OBJECTS, mv changed to /bin/mv
# (to bypass aliasing), cp changed to /bin/cp, add smem and
# testprog targets. See also changes and comments in configure.in


# Default library name:
PACKAGE		= cfitsio

# CFITSIO version numbers:
CFITSIO_MAJOR   = 4
CFITSIO_MINOR	= 1
CFITSIO_MICRO   = 0
CFITSIO_SONAME	= 9

prefix		= /Users/linn/.cfitsio
exec_prefix	= ${prefix}
CFITSIO_BIN	= ${DESTDIR}${exec_prefix}/bin
CFITSIO_LIB	= ${DESTDIR}${exec_prefix}/lib
CFITSIO_INCLUDE	= ${DESTDIR}${prefix}/include
INSTALL_DIRS	= ${DESTDIR}${prefix} ${CFITSIO_INCLUDE} ${CFITSIO_LIB} ${CFITSIO_LIB}/pkgconfig


SHELL =		/bin/sh
ARCHIVE =	ar rv
RANLIB =	ranlib
CC =		gcc
CFLAGS =	-I/Users/linn/miniconda3/include -g -O2 -Dg77Fortran -fPIC -fno-common
CPPFLAGS =	
SSE_FLAGS =	
FC =		gfortran
LDFLAGS =	-g -O2
LDFLAGS_BIN =	-g -O2 -Wl,-rpath,${CFITSIO_LIB}
DEFS =		-DPACKAGE_NAME=\"\" -DPACKAGE_TARNAME=\"\" -DPACKAGE_VERSION=\"\" -DPACKAGE_STRING=\"\" -DPACKAGE_BUGREPORT=\"\" -DPACKAGE_URL=\"\" -DHAVE_STDIO_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_STRINGS_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_UNISTD_H=1 -DSTDC_HEADERS=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MATH_H=1 -DHAVE_LIMITS_H=1 -D_LARGEFILE_SOURCE=1 -D_FILE_OFFSET_BITS=64 -DCFITSIO_HAVE_CURL=1 -DHAVE_LIBZ=1 -DHAVE_FTRUNCATE=1 -DHAVE_LONGLONG=1 -DHAVE_SHMEM_SERVICES=1 -DHAVE_UNION_SEMUN=1 -DHAVE_FMEMOPEN=1 -DHAVE_NET_SERVICES=1
LIBS =		-lz 
LIBS_CURL =	-Wl,-rpath,/Users/linn/miniconda3/lib -L/Users/linn/miniconda3/lib -lcurl
FLEX =		flex
BISON =		bison

SHLIB_LD =	gcc -dynamiclib  -headerpad_max_install_names -install_name @rpath/lib${PACKAGE}.${CFITSIO_SONAME}${SHLIB_SUFFIX} -compatibility_version ${CFITSIO_SONAME} -current_version ${CFITSIO_SONAME}.${CFITSIO_MAJOR}.${CFITSIO_MINOR}.${CFITSIO_MICRO}
SHLIB_SUFFIX =	.dylib
CFITSIO_SHLIB =	lib${PACKAGE}.${CFITSIO_SONAME}.${CFITSIO_MAJOR}.${CFITSIO_MINOR}.${CFITSIO_MICRO}${SHLIB_SUFFIX}
CFITSIO_SHLIB_SONAME = lib${PACKAGE}.${CFITSIO_SONAME}${SHLIB_SUFFIX}


CORE_SOURCES = 	buffers.c cfileio.c checksum.c drvrfile.c drvrmem.c \
		drvrnet.c drvrsmem.c editcol.c edithdu.c eval_l.c \
		eval_y.c eval_f.c fitscore.c getcol.c getcolb.c getcold.c getcole.c \
		getcoli.c getcolj.c getcolk.c getcoll.c getcols.c getcolsb.c \
		getcoluk.c getcolui.c getcoluj.c getkey.c group.c grparser.c \
		histo.c iraffits.c \
		modkey.c putcol.c putcolb.c putcold.c putcole.c putcoli.c \
		putcolj.c putcolk.c putcoluk.c putcoll.c putcols.c putcolsb.c \
		putcolu.c putcolui.c putcoluj.c putkey.c region.c scalnull.c \
		swapproc.c wcssub.c wcsutil.c imcompress.c quantize.c ricecomp.c \
		pliocomp.c fits_hcompress.c fits_hdecompress.c \
		simplerng.c 

ZLIB_SOURCES =	zcompress.c zuncompress.c

SOURCES = ${CORE_SOURCES} ${ZLIB_SOURCES} ${FITSIO_SRC}

OBJECTS = 	${SOURCES:.c=.o}

CORE_OBJECTS = 	${CORE_SOURCES:.c=.o} ${ZLIB_SOURCES:.c=.o}


FITSIO_SRC =	f77_wrap1.c f77_wrap2.c f77_wrap3.c f77_wrap4.c

# ============ description of all targets =============
#       -  <<-- ignore error code

all:
		@if [ "x${FC}" = x ]; then \
			${MAKE} all-nofitsio; \
		else \
			${MAKE} stand_alone; \
		fi

all-nofitsio:
		${MAKE} stand_alone "FITSIO_SRC="

stand_alone:	lib${PACKAGE}.a shared

lib${PACKAGE}.a:	${OBJECTS}
		${ARCHIVE} $@ ${OBJECTS}; \
		${RANLIB} $@;

shared: lib${PACKAGE}${SHLIB_SUFFIX}

lib${PACKAGE}${SHLIB_SUFFIX}: ${OBJECTS}
		${SHLIB_LD} ${LDFLAGS} -o ${CFITSIO_SHLIB} ${OBJECTS} -lm ${LIBS_CURL} ${LIBS}
		@if [ "x${CFITSIO_SHLIB_SONAME}" != x ]; then \
			ln -sf ${CFITSIO_SHLIB} ${CFITSIO_SHLIB_SONAME}; \
			ln -sf ${CFITSIO_SHLIB_SONAME} $@; \
		fi

install:	lib${PACKAGE}.a ${INSTALL_DIRS}
		@for lib in lib${PACKAGE}.a lib${PACKAGE}${SHLIB_SUFFIX} \
				${CFITSIO_SHLIB} ${CFITSIO_SHLIB_SONAME}; do \
		    if [ -f $$lib ]; then \
			echo "/bin/cp -a $$lib ${CFITSIO_LIB}"; \
			/bin/cp -a $$lib ${CFITSIO_LIB}; \
		    fi; \
		done
		/bin/cp fitsio.h fitsio2.h longnam.h drvrsmem.h ${CFITSIO_INCLUDE}
		/bin/cp cfitsio.pc ${CFITSIO_LIB}/pkgconfig
		@for task in ${FPACK_UTILS} ${UTILS}; do \
		    if [ -f $$task ]; then \
			if [ ! -d ${CFITSIO_BIN} ]; then mkdir -p ${CFITSIO_BIN}; fi; \
			echo "/bin/cp $$task ${CFITSIO_BIN}"; \
			/bin/cp $$task ${CFITSIO_BIN}; \
		    fi; \
		done

.c.o:
		${CC} -c -o ${<D}/${@F} ${CFLAGS} ${CPPFLAGS} ${DEFS} $<

UTILS		= cookbook fitscopy imcopy smem speed testprog

FPACK_UTILS	= fpack funpack

# Optional build target:
utils:		${FPACK_UTILS} ${UTILS}

swapproc.o:	swapproc.c
		${CC} -c ${CFLAGS} ${CPPFLAGS} ${SSE_FLAGS} ${DEFS} $<

smem:		smem.o lib${PACKAGE}.a ${OBJECTS}
		${CC} ${LDFLAGS_BIN} ${DEFS} -o $@ ${@}.o -L. -l${PACKAGE} -lm

testprog:	testprog.o lib${PACKAGE}.a ${OBJECTS}
		${CC} ${LDFLAGS_BIN} ${DEFS} -o $@ ${@}.o -L. -l${PACKAGE} -lm ${LIBS}

fpack:		fpack.o fpackutil.o lib${PACKAGE}.a ${OBJECTS}
		${CC} ${LDFLAGS_BIN} ${DEFS} -o $@ ${@}.o fpackutil.o -L. -l${PACKAGE} -lm ${LIBS}

funpack:	funpack.o fpackutil.o lib${PACKAGE}.a ${OBJECTS}
		${CC} ${LDFLAGS_BIN} ${DEFS} -o $@ ${@}.o fpackutil.o -L. -l${PACKAGE} -lm ${LIBS}

fitscopy:	fitscopy.o lib${PACKAGE}.a ${OBJECTS}
		${CC} ${LDFLAGS_BIN} ${DEFS} -o $@ ${@}.o -L. -l${PACKAGE} -lm ${LIBS}

speed:		speed.o lib${PACKAGE}.a ${OBJECTS}
		${CC} ${LDFLAGS_BIN} ${DEFS} -o $@ ${@}.o -L. -l${PACKAGE} -lm ${LIBS}

imcopy:		imcopy.o lib${PACKAGE}.a ${OBJECTS}
		${CC} ${LDFLAGS_BIN} ${DEFS} -o $@ ${@}.o -L. -l${PACKAGE} -lm ${LIBS}

cookbook:	cookbook.o lib${PACKAGE}.a ${OBJECTS}
		${CC} ${LDFLAGS_BIN} ${DEFS} -o $@ ${@}.o -L. -l${PACKAGE} -lm ${LIBS}

testf77:	testf77.o lib${PACKAGE}.a ${OBJECTS}
		${FC} ${LDFLAGS_BIN} -o $@ ${@}.o -L. -l${PACKAGE} -lm ${LIBS}

eval:		# Rebuild eval_* files from flex/bison source
		${FLEX} -t eval.l > eval_l.c1
		# Note workaround for yyfree=fffree conflict
		/bin/sed -e 's/yy/ff/g' -e 's/YY/FF/g' eval_l.c1 -e 's/fffree/yyfffree/g' > eval_l.c
		/bin/rm -f eval_l.c1
		${BISON} -d -v -y eval.y
		/bin/sed -e 's/yy/ff/g' -e 's/YY/FF/g' y.tab.c > eval_y.c
		/bin/sed -e 's/yy/ff/g' -e 's/YY/FF/g' y.tab.h > eval_tab.h
		/bin/rm -f y.tab.c y.tab.h

clean:
	-	/bin/rm -rf *.o zlib/*.o *.dSYM lib${PACKAGE}* y.output so_locations \
		${UTILS} ${FPACK_UTILS} testprog.fit testf77

distclean:	clean
	-	/bin/rm -f Makefile cfitsio.pc config.log config.status configure.lineno

# Make target which outputs the list of the .o contained in the cfitsio lib
# usefull to build a single big shared library containing Tcl/Tk and other
# extensions.  used for the Tcl Plugin. 

cfitsioLibObjs:
	@echo ${CORE_OBJECTS}

cfitsioLibSrcs:
	@echo ${SOURCES}

# This target actually builds the objects needed for the lib in the above
# case
objs: ${CORE_OBJECTS}

${INSTALL_DIRS}:
	@if [ ! -d $@ ]; then mkdir -p $@; fi
