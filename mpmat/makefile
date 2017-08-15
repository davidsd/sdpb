GMPINCLUDEDIR   = /usr/local/include
BOOSTINCLUDEDIR = /usr/local/include/boost
LIBDIR          = /usr/local/lib
MKLROOT         = /opt/intel/mkl

SOURCES = $(wildcard src/*.cpp) $(wildcard src/mpack/*.cpp) $(wildcard src/tests/*.cpp)
HEADERS = $(wildcard src/*.h) $(wildcard src/mpack/*.h) $(wildcard src/tests/*.h)
OBJECTS = $(patsubst src/%.cpp,obj/%.o,$(SOURCES))
RESULT  = mpmat_test

BINARY_DIR = .

CC = /usr/local/opt/llvm/bin/clang


CFLAGS += -I/usr/local/opt/llvm/include \
          -g \
          -Wall -ansi -std=c++0x \
          -Isrc/mpack \
          -I${GMPINCLUDEDIR} \
          -I${BOOSTINCLUDEDIR} \
	  -DMKL_ILP64 -m64 -I${MKLROOT}/include \
          -fopenmp \
          -D___MPACK_BUILD_WITH_GMP___

ifdef MPMAT_DEBUG

    CFLAGS += -O0 -D___MPMAT_DEBUG___

else

    CFLAGS += -O2

endif

LIBS +=  -L/usr/local/opt/llvm/lib \
	 ${MKLROOT}/lib/libmkl_intel_ilp64.a ${MKLROOT}/lib/libmkl_intel_thread.a ${MKLROOT}/lib/libmkl_core.a \
         -liomp5 -lgmpxx -lgmp \
         -lboost_serialization -lboost_system -lboost_filesystem -lboost_timer -lboost_program_options -lboost_chrono \
         -lc++ \
	 -lpthread \
	 -lm \
	 -ldl \

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


