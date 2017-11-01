ifeq ($(HOME), /home/dsd)

GMPINCLUDEDIR   = /usr/include

else

GMPINCLUDEDIR   = /usr/include

endif

all:
	nvcc test2.cu  -I${GMPINCLUDEDIR} -I/group/hepth/rajeeve/bootstrap/boost -lgomp -lgmpxx -lgmp -lcublas -lcurand -lcuda -Xcompiler -fopenmp -DNDEBUG -DBOOST_UBLAS_NDEBUG -G -g -o test
