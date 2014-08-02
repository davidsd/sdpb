OBJECTS = main.o tinyxml2.o
HEADERS = types.h tinyxml2.h
SOURCES = $(OJBECTS:.o=.cpp)
RESULT  = sdp-bootstrap

CC = g++
CFLAGS = -g -O2 -Wall -ansi -L/home/dsd/lib -I/home/dsd/include/mpack -I/home/dsd/include -I/home/dsd/include/boost -fopenmp
RM = rm -f

.SUFFIXES: .cpp .o

$(RESULT): $(OBJECTS)
	$(CC) $(CFLAGS) -lgomp -lmblas_gmp -lmlapack_gmp -lgmp -lgmpxx -lmpc -lboost_serialization -lboost_system -lboost_filesystem -lboost_timer -o $@ $(OBJECTS)

.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	$(RM) *.o *.core core *~ src/*.o math/*.o

$(OBJECTS): $(HEADERS)
