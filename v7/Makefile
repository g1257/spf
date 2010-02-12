# DO NOT EDIT!!! Changes will be lost. Modify configure.pl instead
# This Makefile was written by configure.pl 
# SPF v6.4 by G.A.
# Platform: linux
# RunType for TPEM: serial
# RunType for Trivial: serial
# Model: MODEL_KONDO_PNICTIDES
# GSL: no
LDFLAGS = -L.   -llapack -lblas -lm -L../lib 
EXENAME = spf
CPPFLAGS = -DNDEBUG -I../PartialPsimag   -IGeometries -IModels/PnictidesTwoOrbitals -IEngine -IClassicalFields 
CXX = g++ -Werror -Wall -g3 -pg

$(EXENAME): clean main.o 
	$(CXX) -o $(EXENAME) main.o $(LDFLAGS) 


clean:
	rm -f core* $(EXENAME) *.o *.ii *.tt

$(EXENAME).tar: *.cpp *.h Makefile* *.pl *.inp *.c
	tar -cvf $(EXENAME).tar  *.cpp *.h Makefile* *.pl *.inp *.c

######## End of Makefile ########

