
/** \ingroup SPF */
/*@{*/

/*! \file Spin.h
 *
 *  DynVars: theta and phi (classical spin)
 *
 */
#ifndef SPIN_H
#define SPIN_H
#include "IoSimple.h"

namespace Spf {
	template<typename FieldType_>
	struct Spin { // Do not add functions here, this is a struct!!
		typedef typename PsimagLite::IoSimple::In IoSimpleIn;
		typedef FieldType_ FieldType;
		
		Spin(size_t vol,const std::string& mcstarttype,bool freeze=false,bool makeVoid=false) : 
				size(vol),theta(vol),phi(vol),isFrozen(freeze),isVoid(makeVoid)
		{
			if (mcstarttype=="none") return;
			IoSimpleIn ioin(mcstarttype);
			(*this)<=ioin;
			if (theta.size()==0 || phi.size()==0) throw std::runtime_error("PRoblem\n");
		}
				
		//size_t size() const { return theta.size(); }
		
		size_t size;
		std::vector<FieldType> theta;
		std::vector<FieldType> phi;
		bool isFrozen;
		bool isVoid;
		
	}; // Spin
	
	template<typename FieldType>
	std::ostream& operator<<(std::ostream& os,Spin<FieldType>& dynVars)
	{
		os<<"Theta\n";
		os<<dynVars.theta;
		os<<"Phi\n";
		os<<dynVars.phi;
		os<<"IsFrozen"<<dynVars.isFrozen<<"\n";
		return os;
	}
	
	//! Operator to read Dynvars from file
	template<typename FieldType>
	Spin<FieldType>& operator<=(
			Spin<FieldType>& dynVars,
			typename PsimagLite::IoSimple::In& ioin)
	{
		ioin.read(dynVars.theta,"Theta");
		ioin.read(dynVars.phi,"Phi");
		ioin.readline(dynVars.isFrozen,"IsFrozen");
		
		return dynVars;
	}
	
} // namespace Spf

/*@}*/
#endif
