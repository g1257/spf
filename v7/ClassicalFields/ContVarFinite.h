
/** \ingroup SPF */
/*@{*/

/*! \file ContVarFinite.h
 *
 *  continuous variable finite
 *
 */
#ifndef CONT_VAR_FINITE_H
#define CONT_VAR_FINITE_H
#include "IoSimple.h"
#include "Vector.h"

namespace Spf {
template<typename FieldType_>
struct ContVarFinite { // Do not add functions here, this is a struct!!
	typedef typename PsimagLite::IoSimple::In IoSimpleIn;
	typedef FieldType_ FieldType;
	typedef std::pair<FieldType,FieldType> PairRealType;

	ContVarFinite(size_t vol,const PsimagLite::String& mcstarttype,const PairRealType& bounds1)
	    : size(vol),value(vol, 0.0),bounds(bounds1),isFrozen(false)
	{
		if (mcstarttype=="none") return;
		IoSimpleIn ioin(mcstarttype);

		ioin.read(value,"ContVarFinite");
		ioin.readline(isFrozen,"IsFrozenContVarFinite");

		if (value.size()==0) throw PsimagLite::RuntimeError("Problem in ContVarFinite\n");
	}

	size_t size;
	typename PsimagLite::Vector<FieldType>::Type value;
	PairRealType bounds;
	bool isFrozen;

}; // Spin

template<typename FieldType>
std::ostream& operator<<(std::ostream& os,const ContVarFinite<FieldType>& dynVars)
{
	os<<"ContVarFinite\n";
	os<<dynVars.value;
	os<<"IsFrozenContVarFinite="<<dynVars.isFrozen<<"\n";
	return os;
}

} // namespace Spf

/*@}*/
#endif
