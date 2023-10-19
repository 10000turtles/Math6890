#
#THISMAKEFILEUSESTHEFOLLOWINGENVIRONMENTALVARIABLE:
#APlusPlus
#
#TheAPlusPlusenvvariableshouldnormallybesetwhenyoulogin.
#
#Ifnotyoucansetitusing:
#
#cg6:use:
#setenvAPlusPlus/home/henshw/software/AppPpp-0.8.3-gcc7.5.0/A++/install
#
#cgpi:use:
#setenvAPlusPlus/home/henshw/software/AppPpp-0.8.3-gcc4.8.5/A++/install
#
#Ifyoutype:
#ls$APlusPlus/include 
#youshouldseeA++.haswellasotherinclude files.

#Firsttargetismadebydefaultwhenusing"make",traditionallynamed"all"
all=heatADI

#Setnamesofcompilersonceincaseweneedtochangethem
CC=gcc
CXX=g++

CCFLAGS= -fPIC -O3 -I$(APlusPlus)/include -I../../include

#ListoflibrariesforA++
AppLibraries=-Wl,-rpath,$(APlusPlus)/lib -L$(APlusPlus)/lib -lApp -lApp_static

#Listofalllibraries
LIBS=$(AppLibraries)

#implicit(generic)ruletocompile.Cfiles
#$@=filenameofthetarget
#$<=nameofthefirstprerequistite
%.o : %.C
	$(CXX) $(CCFLAGS) -o $@ -c $<

opt=-O3
FC=gfortran 
FFLAGS= -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -fopenmp
%.o: %.f90 
	$(FC) $(FFLAGS) -o $@ -c $<

#1Dheatequation,implicittime-stepping,A++arrays

heatADIFiles=heatADI.o tridiagonal.o 
heatADI: $(heatADIFiles)
	$(CXX) $(CCFLAGS) -o $@ -fopenmp heatADI.C tridiagonal.o $(LIBS)

heat1dImpFiles=heat2d.o heat2dUpdate.o
heat2d: $(heat1dImpFiles)
	$(CXX) $(CCFLAGS) -o $@ $(heat1dImpFiles) $(LIBS)

poissonFiles = poisson.o
poisson: $(poissonFiles)
	$(CXX) $(CCFLAGS) -o $@ -fopenmp poisson.C $(LIBS)

heatFiles = heat1d.o
heat1d: $(heatFiles)
	$(CXX) $(CCFLAGS) -o  $@ $(heatFiles) $(LIBS)

heat1dp: $(heatFiles)
	$(CXX) $(CCFLAGS) -o $@ -fopenmp heat1d.C  $(LIBS)

heat2dp: $(heat1dImpFiles)
	$(CXX) $(CCFLAGS) -o $@ -fopenmp heat2d.C heat2dUpdate.o $(LIBS)

clean:; rm *.o


