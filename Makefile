UNAME := $(shell uname)

################ Modify these variables ################

ifeq ($(UNAME), Darwin)
# Mac OS X Defaults, assuming Readme instructions.

GMPINCLUDEDIR   = /usr/local/include
BOOSTINCLUDEDIR = /usr/local/include/boost
LIBDIR          = /usr/local/lib

else
# If you're using Linux or Windows or a nonstandard OS X setup, modify
# these variables to point to the appropriate places.

# directory containing the header files 'gmp.h' and 'gmpxx.h'
GMPINCLUDEDIR   = $(HOME)/include

# directory containing boost header files, e.g. 'bind.hpp', etc.
BOOSTINCLUDEDIR = $(HOME)/include/boost

# directory containing library files, e.g.
#   libboost_filesystem.a,
#   libboost_filesystem.so,
#   libboost_filesystem.so.1.54.0
# as well as analogous files for 'boost_system',
# 'boost_serialization', 'boost_timer', 'boost_program_options',
# 'gmp', and 'gmpxx',
LIBDIR          = $(HOME)/lib

endif

################ End of modifications ################

SOURCES := $(wildcard src/*.cpp) $(wildcard src/mpack/*.cpp) $(wildcard src/tinyxml2/*.cpp)
HEADERS := $(wildcard src/*.h) $(wildcard src/mpack/*.h) $(wildcard src/tinyxml2/*.h)
OBJECTS := $(patsubst src/%.cpp,obj/%.o,$(SOURCES))
RESULT  = sdpb

ifdef INTEL

CC = icpc
CFLAGS = -g -O2 -ipo -xhost -Wall -ansi -std=c++0x -L${LIBDIR} -Isrc/mpack -I${GMPINCLUDEDIR} -I${BOOSTINCLUDEDIR} -openmp -D___MPACK_BUILD_WITH_GMP___
LIBS = -lgmpxx -lgmp -lboost_serialization -lboost_system -lboost_filesystem -lboost_timer -lboost_program_options -lboost_chrono -lrt

else 
ifdef CLANG

CC = clang-omp
CFLAGS = -g -O2 -Wall -ansi -std=c++0x -L${LIBDIR} -Isrc/mpack -I${GMPINCLUDEDIR} -I${BOOSTINCLUDEDIR} -fopenmp -D___MPACK_BUILD_WITH_GMP___
LIBS =  -liomp5 -lgmpxx -lgmp -lboost_serialization -lboost_system -lboost_filesystem -lboost_timer -lboost_program_options -lboost_chrono -lc++

else

CC = g++
CFLAGS = -g -O2 -Wall -ansi -std=c++0x -L${LIBDIR} -Isrc/mpack -I${GMPINCLUDEDIR} -I${BOOSTINCLUDEDIR} -fopenmp -D___MPACK_BUILD_WITH_GMP___
LIBS = -lgomp -lgmpxx -lgmp -lboost_serialization -lboost_system -lboost_filesystem -lboost_timer -lboost_program_options -lboost_chrono -lrt

endif
endif


.SUFFIXES: .cpp .o

$(RESULT): $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

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

$(OBJECTS): $(HEADERS) | obj

CFLAGS += -MMD
-include $(OBJECTS:.o=.d)
