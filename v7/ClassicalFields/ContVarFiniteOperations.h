/** \ingroup SPF */
/*@{*/

/*!
 *
 *
 *
 */

#ifndef CONT_VAR_FINITE_OPS_H
#define CONT_VAR_FINITE_OPS_H

#include "Vector.h"
#include "ContVarFinite.h"
#include "ProgramGlobals.h"

namespace Spf {
template<typename GeometryType,typename FieldType>
class ContVarFiniteOperations {

	typedef ContVarFinite<FieldType> ContVarFiniteType;
	typedef typename ContVarFiniteType::PairRealType PairRealType;

public:

	typedef ContVarFiniteType DynVarsType;

	ContVarFiniteOperations(const GeometryType& geometry,
	                        const FieldType& mcwindow,
	                        const PairRealType& bounds)
	    : geometry_(geometry),
	      mcwindow_(mcwindow),
	      dynVars2_(0,"none",bounds)
	{
	}

	void set(DynVarsType& dynVars)
	{
		dynVars_=&dynVars;
	}

	//! How to sweep the lattice
	template<typename RngType>
	size_t proposeSite(size_t i,RngType& rng) const
	{
		return i; //<-- zig-zag horizontal
		// zig-zag vertical:
		/*size_t l = geometry_.length();
			size_t x = i % l;
			size_t y = i / l;
			return y + x*l;*/
		// random:
		//return size_t(rng()*geometry_.volume());

	}

	template<typename RngType>
	void proposeChange(size_t i,RngType& rng)
	{
		FieldType phononsOld = dynVars_->value[i];

		dynVars2_ = *dynVars_;

		propose_(phononsOld,dynVars2_.value[i],rng);
	}

	const DynVarsType& dynVars2() const { return dynVars2_; }

	FieldType deltaDirect(size_t i,const FieldType& coupling) const
	{
		return dSDirect(*dynVars_,dynVars2_,i,coupling);
	}

	FieldType sineUpdate(size_t i) const
	{
		return 1.0; // measure
	}

	void accept(size_t i)
	{
		dynVars_->value[i]=dynVars2_.value[i];
	}
	
	FieldType sum2() const
	{
		FieldType sum = 0;
		for (size_t i=0;i<dynVars_->value.size();i++)
			sum += ProgramGlobals::square(dynVars_->value[i]);
		return sum;
	}


	FieldType sum() const
	{
		FieldType sum = 0;
		for (size_t i=0;i<dynVars_->value.size();i++)
			sum += dynVars_->value[i];
		return sum;
	}

	const PairRealType& bounds() const
	{
		return dynVars2_.bounds;
	}

private:

	const GeometryType& geometry_;
	const FieldType& mcwindow_;
	DynVarsType* dynVars_;
	DynVarsType dynVars2_;

	template<typename RngType>
	void propose_(
	        const FieldType& phononsOld,
	        FieldType& phononsNew,
	        RngType& rng)
	{

		phononsNew=phononsOld + (rng.random()- 0.5)*mcwindow_;
		FieldType lowBound = bounds().first;
		FieldType highBound = bounds().second;
		if (phononsNew>highBound) phononsNew = highBound;
		if (phononsNew<lowBound) phononsNew = lowBound;
	}

	FieldType dSDirect(const DynVarsType& dynVars,const DynVarsType& dynVars2, size_t i,
	                   const FieldType& coupling) const
	{
		FieldType tmp = ProgramGlobals::square(dynVars2.value[i]);
		tmp -= ProgramGlobals::square(dynVars.value[i]);
		for (size_t k=0;k<geometry_.z(1);k++) {
			size_t j = geometry_.neighbor(i,k).first;
			tmp += ProgramGlobals::square(dynVars2.value[j]);
			tmp -= ProgramGlobals::square(dynVars.value[j]);
		}
		return coupling*tmp;
	}
}; // ContVarFiniteOperations

} // namespace Spf

/*@}*/
#endif
