### MKL ###
# Information on how to link MKL can be found on 
# https://software.intel.com/content/www/us/en/develop/articles/intel-mkl-link-line-advisor.html

# Debian with package intel-mkl
MKLLIBDIR=/usr/lib/x86_64-linux-gnu/
MKL_CFLAGS=-DMKL_ILP64 -m64 -I/usr/include/mkl
MKL_LDFLAGS=-Wl,--start-group $(MKLLIBDIR)/libmkl_intel_ilp64.a $(MKLLIBDIR)/libmkl_sequential.a $(MKLLIBDIR)/libmkl_core.a -Wl,--end-group -lm -ldl

# Custom installation in /opt/aponcet/intel/ (on gambit)
#MKLROOT=/opt/aponcet/intel/system_studio_2018/compilers_and_libraries_2018.2.199/linux/mkl/
# Or custom installation in /opt/intel-18/ (on ponyo or chihiro)
#MKLROOT=/opt/intel-18/compilers_and_libraries/linux/mkl/
#MKL_CFLAGS=-DMKL_ILP64 -m64 -I${MKLROOT}/include
#MKL_LDFLAGS=-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group
### END MKL ###

CC=g++
CFLAGS=-W -Wall -std=c++11 -Ofast -march=native $(MKL_CFLAGS) 
#CC=icc
#CFLAGS=-W -Wall -xSSE2 -axAVX,CORE-AVX2,CORE-AVX512 -Ofast $(MKL_CFLAGS) $(MPI_CFLAGS)
LDFLAGS=-lboost_program_options -lpthread $(MKL_LDFLAGS)
EXEC=Comb
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
