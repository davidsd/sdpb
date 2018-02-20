UNAME := $(shell uname)

################ Modify these variables ################

ifeq ($(UNAME), Darwin)
# Mac OS X Defaults, assuming Readme instructions.

GMPINCLUDEDIR   = /usr/local/include
BOOSTINCLUDEDIR = /group/hepth/rajeeve/bootstrap/boost/
LIBDIR          = /usr/local/lib/group/hepth/rajeeve/bootstrap/boost/stage/lib/

else
# If you're using Linux or Windows or a nonstandard OS X setup, modify
# these variables to point to the appropriate places.

BASEDIR := $(HOME)

# directory containing the header files 'gmp.h' and 'gmpxx.h'
GMPINCLUDEDIR   = /usr/include

# directory containing boost header files, e.g. 'bind.hpp', etc.
BOOSTINCLUDEDIR = /home/rajeeve/boost/include
BOOSTLIBDIR = /home/rajeeve/boost/lib
# directory containing library files, e.g.
#   libboost_filesystem.a,
#   libboost_filesystem.so,
#   libboost_filesystem.so.1.54.0
# as well as analogous files for 'boost_system',
# 'boost_serialization', 'boost_timer', 'boost_program_options',
# 'gmp', and 'gmpxx',
LIBDIR          = /usr/lib64
#INTEL 		= yes
endif

MKLROOT         = /opt/intel/compilers_and_libraries_2017.4.196/linux/mkl

CUDAROOT	= /usr/local/cuda-7.5/lib64

################ End of modifications ################

SOURCES = $(wildcard src/*.cpp) $(wildcard src/mpack/*.cpp) $(wildcard src/mpmat/*.cpp)
CUDASOURCES = $(wildcard src/*.cu)
HEADERS = $(wildcard src/*.h) $(wildcard src/mpack/*.h) $(wildcard src/mpmat/*.h)
OBJECTS = $(patsubst src/%.cpp,obj/%.o,$(SOURCES)) 
CUDAOBJECTS = $(patsubst src/%.cu,obj/%.co,$(CUDASOURCES))
RESULT  = sdpb

ifeq ($(SHARED_TINYXML2), 1)
LIBS = -ltinyxml2
CFLAGS += -D___SHARED_TINYXML2___
else
SOURCES += $(wildcard src/tinyxml2/*.cpp)
HEADERS += $(wildcard src/tinyxml2/*.h)
endif

ifdef CLION_EXE_DIR
	BINARY_DIR=$(CLION_EXE_DIR)
	CLANG=1
else
	BINARY_DIR=.
endif

ifdef INTEL

CC = icpc
CFLAGS += -O3 -ipo -xhost -Wall -ansi -mkl -std=c++14 -L${LIBDIR} -DMKL_ILP64 -m64 -I${MKLROOT}/include -Isrc/mpack -I${BOOSTINCLUDEDIR} -I${GMPINCLUDEDIR} -qopenmp -D___MPACK_BUILD_WITH_GMP___
LIBS += ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.so ${MKLROOT}/lib/intel64/libmkl_gnu_thread.so ${MKLROOT}/lib/intel64/libmkl_core.so -lgmpxx -lgmp ${BOOSTLIBDIR}/libboost_serialization.so ${BOOSTLIBDIR}/libboost_system.so ${BOOSTLIBDIR}/libboost_filesystem.so ${BOOSTLIBDIR}/libboost_timer.so ${BOOSTLIBDIR}/libboost_program_options.so ${BOOSTLIBDIR}/libboost_chrono.so -lrt

else 
ifdef CLANG

CC = /usr/local/opt/llvm/bin/clang
CFLAGS += -I/usr/local/opt/llvm/include -O3 -Wall -ansi -std=c++0x -Isrc/mpack -I${GMPINCLUDEDIR} -I${BOOSTINCLUDEDIR} -fopenmp -D___MPACK_BUILD_WITH_GMP___
LIBS +=  -L/usr/local/opt/llvm/lib -liomp5 -lgmpxx -lgmp -lboost_serialization -lboost_system -lboost_filesystem -lboost_timer -lboost_program_options -lboost_chrono -lc++

else

CC = g++
CFLAGS += -O3 -Wall -ansi -std=c++0x -pg -L${LIBDIR} -DMKL_ILP64 -m64 -I${MKLROOT}/include  -Isrc/mpack -I${GMPINCLUDEDIR} -I${BOOSTINCLUDEDIR} -fopenmp -D___MPACK_BUILD_WITH_GMP___
LIBS += ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.so ${MKLROOT}/lib/intel64/libmkl_gnu_thread.so ${MKLROOT}/lib/intel64/libmkl_core.so -lgomp -lgmpxx -lgmp ${BOOSTLIBDIR}/libboost_serialization.so ${BOOSTLIBDIR}/libboost_system.so ${BOOSTLIBDIR}/libboost_filesystem.so ${BOOSTLIBDIR}/libboost_timer.so ${BOOSTLIBDIR}/libboost_program_options.so ${BOOSTLIBDIR}/libboost_chrono.so
#LIBS += ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.so ${MKLROOT}/lib/intel64/libmkl_gnu_thread.so ${MKLROOT}/lib/intel64/libmkl_core.so -lgomp -lgmpxx -lgmp -lboost_serialization -lboost_system -lboost_filesystem -lboost_timer -lboost_program_options -lboost_chrono -lrt 

endif
endif


cuda: LIBS += -lcublas -lcurand -lcudart
cuda: CFLAGS += -D__SDPB_CUDA__

.PHONY: cuda

.SUFFIXES: .cpp .o .cu

$(RESULT): $(OBJECTS) $(CUDAOBJECTS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $(BINARY_DIR)/$@ $^ $(LIBS)

obj/%.o: src/%.cpp
	$(CC) $(CFLAGS) -c -o $@ $<

obj/%.co: src/%.cu
	nvcc $< -c -I${GMPINCLUDEDIR} -I/group/hepth/rajeeve/bootstrap/boost -lgomp -lgmpxx -lgmp -lcublas -lcurand -lcuda -Xcompiler -fopenmp -DNDEBUG -DBOOST_UBLAS_NDEBUG -G -g -o $@

clean:
	rm -rf obj

cuda: $(OBJECTS) $(CUDAOBJECTS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $(BINARY_DIR)/$(RESULT) $^ $(LIBS)

obj:
	@mkdir -p $@
	@mkdir -p $@/mpack
	@mkdir -p $@/tinyxml2
	@mkdir -p $@/mpmat

test:
	./sdpb -s test.xml --noFinalCheckpoint --maxThreads=16 --multTest
	test -f test.out 
	grep "primalObjective = 1.84026576313204924668804017173" test.out > /dev/null || (echo "Test failed: primalObjective differs from 1.84026576313204924668804017173." && exit 1)

$(OBJECTS): $(HEADERS) | obj

CFLAGS += -MMD
-include $(OBJECTS:.o=.d)
