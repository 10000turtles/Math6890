# export APlusPlus=/home/henshw/software/AppPpp-0.8.3-gcc7.5.0-parallel/A++/install
# export mpi=/home/henshw/software/mpich-3.3.1-install

APlusPlus = /gpfs/u/barn/PCM1/shared/software/AppPpp-0.8.3-parallel/P++/install

all=openMPIHW

#Setnamesofcompilersonceincaseweneedtochangethem
CC=gcc
CXX=g++
CMPI=/home/henshw/software/mpich-3.3.1-install/bin/mpic++

CCFLAGS= -DUSE_PPP -fPIC -O3 -I$(mpi)/include -I$(APlusPlus)/include -I../HW1thru/include

#ListoflibrariesforA++
AppLibraries=-Wl,-rpath,$(APlusPlus)/lib -L$(APlusPlus)/lib -lPpp -lPpp_static

LIBS=$(AppLibraries)

%.o : %.C
	$(CXX) $(CCFLAGS) -o $@ -c $<

opt=-O3
FC=gfortran 
FFLAGS= -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -fopenmp
%.o: %.f90 
	$(FC) $(FFLAGS) -o $@ -c $<

openMPIFiles = heat2d.o heat2dUpdate.o
heat2d: $(openMPIFiles)
	$(CMPI) $(CCFLAGS) -o $@ $(openMPIFiles) $(LIBS)


clean:; rm *.o
