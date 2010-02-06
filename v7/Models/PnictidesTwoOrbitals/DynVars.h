
/** \ingroup SPF */
/*@{*/

/*! \file DynVars.h
 *
 *  DynVars: theta and phi (classical fields)
 *
 */
#ifndef DYNVARS_H
#define DYNVARS_H
#include "Utils.h"

namespace Spf {
	template<typename FieldType>
	struct DynVars { // Do not add functions here, this is a struct!!
		std::vector<FieldType> theta;
		std::vector<FieldType> phi;
		bool isFrozen;
	}; // DynVars
} // namespace Spf

/*@}*/
#endif
