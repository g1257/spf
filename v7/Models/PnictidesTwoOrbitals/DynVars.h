
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
		DynVars(size_t vol) : theta(vol),phi(vol),isFrozen(false)
		{} 
				
		std::vector<FieldType> theta;
		std::vector<FieldType> phi;
		bool isFrozen;
	}; // DynVars
	
	template<typename FieldType>
	std::ostream& operator<<(std::ostream& os,DynVars<FieldType>& dynVars)
	{
		os<<"Theta\n";
		os<<dynVars.theta;
		os<<"Phi\n";
		os<<dynVars.phi;
		os<<"IsFrozen "<<dynVars.isFrozen<<"\n";
		return os;
	}
	
} // namespace Spf

/*@}*/
#endif
