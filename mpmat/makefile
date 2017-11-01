GMPINCLUDEDIR   = /usr/local/include
BOOSTINCLUDEDIR = /group/hepth/rajeeve/bootstrap/boost/
#BOOSTINCLUDEDIR = /home/liliesiu/numerics/boost/
#LIBDIR          = /usr/local/lib
LIBDIR          = /group/hepth/rajeeve/bootstrap/boost/stage/lib/
MKLROOT         = /opt/intel/compilers_and_libraries_2017.4.196/linux/mkl

SOURCES = $(wildcard src/*.cpp) $(wildcard src/mpack/*.cpp) $(wildcard src/tests/*.cpp)
HEADERS = $(wildcard src/*.h) $(wildcard src/mpack/*.h) $(wildcard src/tests/*.h)
OBJECTS = $(patsubst src/%.cpp,obj/%.o,$(SOURCES))
RESULT  = mpmat_test

MPMAT_DEBUG = ""

BINARY_DIR = .

CC = g++


CFLAGS += -g \
          -Wall -ansi -std=c++0x\
          -Isrc/mpack \
          -I${BOOSTINCLUDEDIR} \
	  -I${GMPINCLUDEDIR} \
	  -L${LIBDIR}\
	  -DMKL_ILP64 -m64 -I${MKLROOT}/include \
          -fopenmp \
          -D___MPACK_BUILD_WITH_GMP___

ifdef MPMAT_DEBUG

    CFLAGS += -O0 -D___MPMAT_DEBUG___

else

    CFLAGS += -O2

endif

#${MKLROOT}/lib/intel64/libmkl_gnu_thread.a

LIBS +=	 ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.so ${MKLROOT}/lib/intel64/libmkl_gnu_thread.so ${MKLROOT}/lib/intel64/libmkl_core.so \
         -lgomp -lgmpxx -lgmp \
	 -lpthread \
	 -lm \
	 -ldl \
        -lboost_serialization -lboost_system -lboost_filesystem -lboost_timer -lboost_program_options -lboost_chrono \
	#/group/hepth/rajeeve/bootstrap/boost/stage/lib/libboost_serialization.so.1.61.0 -lboost_system -lboost_filesystem -lboost_timer -lboost_program_options -lboost_chrono \
         #


.SUFFIXES: .cpp .o

$(RESULT): $(OBJECTS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $(BINARY_DIR)/$@ $^ $(LIBS)
	cp $@ cmake-build-debug/mpmat

run : $(RESULT)
	./$(RESULT)

obj/%.o: src/%.cpp
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -rf obj

obj:
	@mkdir -p $@
	@mkdir -p $@/mpack
	@mkdir -p $@/tinyxml2
	@mkdir -p $@/tests

$(OBJECTS): $(HEADERS) | obj

CFLAGS += -MMD
-include $(OBJECTS:.o=.d)


