%module myspf
%{
	#include "GeometrySquare.h"
	#include "ConcurrencySerial.h"
%}

namespace Spf {

template<typename T>
class GeometrySquare {
	public:
	GeometrySquare(size_t l);
	size_t length() const;
	std::string name() const;
};
%template(GeometrySquareD) GeometrySquare<double>;
} // namespace Spf 

namespace Dmrg {
template<typename FieldType>
class ConcurrencySerial : public Concurrency<FieldType> {
	public:
		
	ConcurrencySerial(int argc,char *argv[]);
		
	int nprocs();
		
	int rank();
};
%template(ConcurrencySerialD) ConcurrencySerial<double>;	


} // namespace Dmrg