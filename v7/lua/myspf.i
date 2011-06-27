%module myspf
%{
	#include "GeometrySquare.h"
%}

namespace Spf {

template<typename T>
class GeometrySquare {
	public:
	GeometrySquare(size_t l);
};


%template(GeometrySquareD) GeometrySquare<double>;
} // namespace Spf