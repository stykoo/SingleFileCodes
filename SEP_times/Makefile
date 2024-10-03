### MKL ###
# Information on how to link MKL can be found on 
# https://software.intel.com/content/www/us/en/develop/articles/intel-mkl-link-line-advisor.html

MKLROOT=/opt/intel/mkl/
MKL_CFLAGS=-DMKL_ILP64 -m64 -I${MKLROOT}/include
MKL_LDFLAGS=${MKLROOT}/lib/libmkl_intel_ilp64.a ${MKLROOT}/lib/libmkl_sequential.a ${MKLROOT}/lib/libmkl_core.a -lm -ldl
### END MKL ###

CC=g++
CFLAGS=-g -W -Wall -Ofast -std=c++11 -march=native $(MKL_CFLAGS) 
#CC=icc
#CFLAGS=-W -Wall -xSSE2 -axAVX,CORE-AVX2,CORE-AVX512 -Ofast $(MKL_CFLAGS) $(MPI_CFLAGS)
LDFLAGS=-g -lboost_program_options -lpthread $(MKL_LDFLAGS)
EXEC=SEP_times
SRC=$(wildcard *.cpp)
OBJ=$(SRC:.cpp=.o)

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

main.o: main.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

%.o: %.cpp %.h parameters.h
	$(CC) -o $@ -c $< $(CFLAGS)

.PHONY: clean mrproper

clean:
		rm -f *.o

mrproper: clean
		rm -rf $(EXEC)
