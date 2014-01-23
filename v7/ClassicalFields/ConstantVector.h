
/** \ingroup SPF */
/*@{*/

/*! \file ConstantVector.h
 *
 *
 *
 */
#ifndef CONSTANT_VECTOR_H
#define CONSTANT_VECTOR_H

#include "Vector.h"

namespace Spf {

class ConstantVector {

public:

	SizeType operator[](SizeType i) { return 1; }
};

} // namespace Spf

/*@}*/
#endif
