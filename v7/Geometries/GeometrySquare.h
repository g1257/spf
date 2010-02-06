
/** \ingroup SPF */
/*@{*/

/*! \file GeometrySquare.h
 *
 * A square lattice
 *
 */
#ifndef GEOM_SQUARE_H
#define GEOM_SQUARE_H
#include "Utils.h"

namespace Spf {
	template<typename FieldType_>
	class GeometrySquare {
		public:
		//typedef FieldType_ FieldType;
		typedef std::pair<size_t,size_t> PairType;
		
		GeometrySquare(size_t l) : l_(l) 
		{
			volume_ = l*l;
		}
		
		size_t volume() const { return volume_; }
		
		private:
		
		size_t l_;
		size_t volume_;
	};
	
} // namespace Spf


/*@}*/
#endif
