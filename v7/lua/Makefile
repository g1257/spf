
CXXFLAGS= -I../Geometries -I../Models/PnictidesMultiOrbitals -I../Engine -I../ClassicalFields -I../../../PsimagLite/src

all: myspf.so


myspf.so: myspf_wrap.o  
	g++ $(CXXFLAGS)  -shared myspf_wrap.o -o myspf.so -llua -llapack -lblas

myspf_wrap.o: myspf_wrap.cxx
	g++ $(CXXFLAGS) -static -fPIC -c myspf_wrap.cxx 

myspf_wrap.cxx: myspf.i 
	swig -c++ -lua $(CXXFLAGS) myspf.i

clean:
	rm -f *.o min.o *.so

 
