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
GMPINCLUDEDIR   = $(BASEDIR)/include

# directory containing boost header files, e.g. 'bind.hpp', etc.
BOOSTINCLUDEDIR =  /group/hepth/rajeeve/bootstrap/boost/

# directory containing library files, e.g.
#   libboost_filesystem.a,
#   libboost_filesystem.so,
#   libboost_filesystem.so.1.54.0
# as well as analogous files for 'boost_system',
# 'boost_serialization', 'boost_timer', 'boost_program_options',
# 'gmp', and 'gmpxx',
LIBDIR          = /group/hepth/rajeeve/bootstrap/boost/stage/lib/

endif

MKLROOT         = /opt/intel/compilers_and_libraries_2017.4.196/linux/mkl

################ End of modifications ################

SOURCES = $(wildcard src/*.cpp) $(wildcard src/mpack/*.cpp)
HEADERS = $(wildcard src/*.h) $(wildcard src/mpack/*.h)
OBJECTS = $(patsubst src/%.cpp,obj/%.o,$(SOURCES))
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
CFLAGS += -g -O0 -ipo -xhost -Wall -ansi -std=c++0x -L${LIBDIR} -Isrc/mpack -I${GMPINCLUDEDIR} -I${BOOSTINCLUDEDIR} -openmp -D___MPACK_BUILD_WITH_GMP___
LIBS += -lgmpxx -lgmp -lboost_serialization -lboost_system -lboost_filesystem -lboost_timer -lboost_program_options -lboost_chrono -lrt

else 
ifdef CLANG

CC = /usr/local/opt/llvm/bin/clang
CFLAGS += -I/usr/local/opt/llvm/include -g -O0 -Wall -ansi -std=c++0x -Isrc/mpack -I${GMPINCLUDEDIR} -I${BOOSTINCLUDEDIR} -fopenmp -D___MPACK_BUILD_WITH_GMP___
LIBS +=  -L/usr/local/opt/llvm/lib -liomp5 -lgmpxx -lgmp -lboost_serialization -lboost_system -lboost_filesystem -lboost_timer -lboost_program_options -lboost_chrono -lc++

else

CC = g++
CFLAGS += -O3 -Wall -ansi -std=c++0x -L${LIBDIR} -DMKL_ILP64 -m64 -I${MKLROOT}/include  -Isrc/mpack -I${GMPINCLUDEDIR} -I${BOOSTINCLUDEDIR} -fopenmp -D___MPACK_BUILD_WITH_GMP___
LIBS += ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.so ${MKLROOT}/lib/intel64/libmkl_gnu_thread.so ${MKLROOT}/lib/intel64/libmkl_core.so -lgomp -lgmpxx -lgmp -lboost_serialization -lboost_system -lboost_filesystem -lboost_timer -lboost_program_options -lboost_chrono -lrt

endif
endif


.SUFFIXES: .cpp .o


$(RESULT): $(OBJECTS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $(BINARY_DIR)/$@ $^ $(LIBS)

obj/%.o: src/%.cpp
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -rf obj

obj:
	@mkdir -p $@
	@mkdir -p $@/mpack
	@mkdir -p $@/tinyxml2

test:
	./sdpb -s test.xml --noFinalCheckpoint
	test -f test.out
	grep "primalObjective = 1.84026576313204924668804017173" test.out > /dev/null || \
	(echo "Test failed: primalObjective differs from 1.84026576313204924668804017173." && exit 1)

$(OBJECTS): $(HEADERS) | obj

CFLAGS += -MMD
-include $(OBJECTS:.o=.d)
