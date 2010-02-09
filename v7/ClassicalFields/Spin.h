
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
#include "IoSimple.h"

namespace Spf {
	template<typename FieldType_>
	struct DynVars { // Do not add functions here, this is a struct!!
		typedef typename Dmrg::IoSimple::In IoSimpleIn;
		typedef FieldType_ FieldType;
		
		DynVars(size_t vol,const std::string& mcstarttype) : theta(vol),phi(vol),isFrozen(false)
		{
			if (mcstarttype=="none") return;
			IoSimpleIn ioin(mcstarttype);
			(*this)<=ioin;
			if (theta.size()==0 || phi.size()==0) throw std::runtime_error("PRoblem\n");
		}
				
		size_t size() const { return theta.size(); }
		
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
	
	//! Operator to read Dynvars from file
	template<typename FieldType>
	DynVars<FieldType>&
	operator <= (DynVars<FieldType>& dynVars,  typename Dmrg::IoSimple::In& ioin) 
	{
		ioin.read(dynVars.theta,"Theta");
		ioin.read(dynVars.phi,"Phi");
		ioin.readline(dynVars.isFrozen,"IsFrozen");
		
		return dynVars;
	}
	
} // namespace Spf

/*@}*/
#endif
