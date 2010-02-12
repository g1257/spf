# DO NOT EDIT!!! Changes will be lost. Modify driver.pl instead
# This Makefile was written by driver.pl
# SPF v7 by G.A. and C.S.
# Platform: linux
# Model: PhononsTwoOrbitals
LDFLAGS = -L.   -llapack -lblas -lm -L../lib 
EXENAME = spf
CPPFLAGS = -DNDEBUG -I../PartialPsimag   -IGeometries -IModels/PhononsTwoOrbitals -IEngine -IClassicalFields 
CXX = g++ -Werror -Wall -g3 -pg

$(EXENAME): clean main.o 
	$(CXX) -o $(EXENAME) main.o $(LDFLAGS) 


clean:
	rm -f core* $(EXENAME) *.o *.ii *.tt

######## End of Makefile ########

