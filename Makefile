OBJECTS = main.o mpack/Rdot.o \
        mpack/Rtrsm.o mpack/Rsyrk.o mpack/Raxpy.o \
        mpack/Rgemm.o mpack/Rtrmm.o mpack/Rtrsv.o \
        mpack/iMlaenv.o mpack/Rlamch.o mpack/Rlascl.o \
        mpack/Rsytrd.o mpack/Rsterf.o mpack/Rorgtr.o \
        mpack/Rlatrd.o mpack/Rsyr2k.o mpack/Rsytd2.o \
        mpack/Rlanst.o mpack/Rlae2.o mpack/Rlapy2.o \
        mpack/Rlasrt.o mpack/Rorgql.o mpack/Rorgqr.o \
        mpack/Rsymv.o mpack/Rlarfg.o mpack/Rsyr2.o \
        mpack/Rlassq.o mpack/Rorg2l.o mpack/Rlarft.o \
        mpack/Rlarfb.o mpack/Rorg2r.o mpack/Rnrm2.o \
        mpack/Rlarf.o mpack/Rger.o mpack/Rpotrf.o \
        mpack/Mxerbla.o mpack/Rpotf2.o mpack/Mlsame.o \
        mpack/Rscal.o mpack/Rcopy.o mpack/Rgemv.o \
        mpack/Rtrmv.o mpack/Rsteqr.o mpack/Rlaset.o \
        mpack/Rlaev2.o mpack/Rlasr.o mpack/Rlartg.o \
        mpack/Rswap.o mpack/Rsyev.o mpack/Rlansy.o \
        mpack/Mutils.o
HEADERS = types.h mpack/mblas_gmp.h mpack/mlapack_gmp.h mpack/mpack_config.h mpack/mutils_gmp.h SquareMatrix.h
SOURCES = $(OJBECTS:.o=.cpp)
RESULT  = sdp-bootstrap

CC = g++
CFLAGS = -g -O2 -Wall -ansi -pedantic -L/home/dsd/lib -I./mpack -I/home/dsd/include -I/home/dsd/include/boost -fopenmp
RM = rm -f

.SUFFIXES: .cpp .o

$(RESULT): $(OBJECTS)
	$(CC) $(CFLAGS) -lgomp -lmblas_gmp -lmlapack_gmp -lgmp -lgmpxx -lmpc -lboost_serialization -lboost_system -lboost_filesystem -o $@ $(OBJECTS)

.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	$(RM) *.o *.core core *~ src/*.o math/*.o

$(OBJECTS): $(HEADERS)
