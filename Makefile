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

CCFLAGS= -fPIC -O1 -I$(APlusPlus)/include -I../../include

#ListoflibrariesforA++
AppLibraries=-Wl,-rpath,$(APlusPlus)/lib -L$(APlusPlus)/lib -lApp -lApp_static

#Listofalllibraries
LIBS=$(AppLibraries)

#implicit(generic)ruletocompile.Cfiles
#$@=filenameofthetarget
#$<=nameofthefirstprerequistite
%.o : %.C
	$(CXX) $(CCFLAGS) -o $@ -c $<

opt=-01
FC=gfortran 
FFLAGS=$(opt)-fdefault-real-8 -fdefault-double-8 -ffree-line-length-none 
%.o: %.f90 
	$(FC) $(FFLAGS) -o $@ -c $<

#1Dheatequation,implicittime-stepping,A++arrays

heatADIFiles=heatADI.o tridiagonal.o 
heatADI: $(heatADIFiles)
	$(CXX) $(CCFLAGS) -o $@ $(heatADIFiles) $(LIBS)

heat1dImpFiles=heat2d.o heat2dUpdate.o
heat2d: $(heat1dImpFiles)
	$(CXX) $(CCFLAGS) -o $@ $(heat1dImpFiles) $(LIBS)





clean:; rm *.o


