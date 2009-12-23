# DO NOT EDIT!!! Changes will be lost. Modify configure.pl instead
# This Makefile was written by configure.pl 
# SPF v6.4 by G.A.
# Platform: linux
# RunType for TPEM: serial
# RunType for Trivial: serial
# Model: MODEL_KONDO_INF_TWOBANDS
# GSL: no
LDFLAGS = -L.   -llapack -lblas -lm  -L../lib 
EXENAME = spf2b
CPPFLAGS = -I../PartialPsimag    -DMODEL_KONDO_INF_TWOBANDS -DNO_GSL  -Wall -DNDEBUG 
CXX = g++ -g3 
include filelist.make

OBJECTS += kondo_inf_twobands.o
OBJECTS += serial_launcher.o
all: $(EXENAME)

tpemSample: tpemSample.o tpemplus.o tpemplusSparse.o tpemplusSubspace.o tpemplus.h tpemplusTypes.h
	$(CXX) tpemSample.o tpemplus.o tpemplusSparse.o tpemplusSubspace.o -o tpemSample $(LDFLAGS)


$(EXENAME): $(HEADERS) $(OBJECTS)
	$(CXX) -o $(EXENAME) $(OBJECTS) $(LDFLAGS)


clean:
	rm -f core* $(EXENAME) *.o *.ii *.tt

$(EXENAME).tar: *.cpp *.h Makefile* *.pl *.inp *.c
	tar -cvf $(EXENAME).tar  *.cpp *.h Makefile* *.pl *.inp *.c

######## End of Makefile ########

