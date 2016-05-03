GMPINCLUDEDIR   = /usr/local/include

all:
	nvcc test.cu -O3 -I${GMPINCLUDEDIR} -lgomp -lgmpxx -lgmp -lcublas -lcurand -lcuda -Xcompiler -fopenmp -DNDEBUG -DBOOST_UBLAS_NDEBUG -o test