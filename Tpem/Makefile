# TPEM by G.A. 
# Platform: linux
LDFLAGS =  -llapack -lblas -lm  -lgsl -lgslcblas 
EXENAME = tpemSample 
CPPFLAGS =  -I../../PsimagLite/src    
CXX = g++ -DNDEBUG -Werror -Wall -O3 -DUSE_GSL
#comment out this one for debugging
#CXX = g++ -Werror -Wall -g3 -DUSE_GSL

all: $(EXENAME)

tpemSample: tpemSample.o
	$(CXX) -o tpemSample tpemSample.o $(LDFLAGS) 

# dependencies brought about by Makefile.dep
%.o: %.cpp Makefile
	$(CXX) $(CPPFLAGS) -c $<

Makefile.dep: tpemSample.cpp Makefile
	$(CXX) $(CPPFLAGS) -MM tpemSample.cpp  > Makefile.dep

clean:
	rm -f core* $(EXENAME) *.o *.ii *.tt

include Makefile.dep
######## End of Makefile ########

