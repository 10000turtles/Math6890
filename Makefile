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
all=heat1dImp

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

#1Dheatequation,implicittime-stepping,A++arrays
heat1dImpFiles=heat1dImp.o tridiagonal.o
heat1dImp: $(heat1dImpFiles)
	$(CXX) $(CCFLAGS) -o $@ $(heat1dImpFiles) $(LIBS)

run:
	./heat1dImp -sol=true -Nx=160 -bc1=n -bc2=n -matlabFileName=heat1d10.m
	./heat1dImp -sol=true -Nx=160 -bc1=n -bc2=n -matlabFileName=heat1d20.m
	./heat1dImp -sol=true -Nx=160 -bc1=n -bc2=n -matlabFileName=heat1d40.m
	./heat1dImp -sol=true -Nx=160 -bc1=n -bc2=n -matlabFileName=heat1d80.m
	./heat1dImp -sol=true -Nx=160 -bc1=n -bc2=n -matlabFileName=heat1d160.m

clean:; rm *.o