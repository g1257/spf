
/** \ingroup SPF */
/*@{*/

/*! \file Phonon.h
 *
 *  DynVars: classical phonons
 *
 */
#ifndef PHONON_H
#define PHONON_H
#include "Utils.h"
#include "IoSimple.h"
#include "Vector.h"

namespace Spf {
	template<typename FieldType_>
	struct Phonon { // Do not add functions here, this is a struct!!
		typedef typename Dmrg::IoSimple::In IoSimpleIn;
		typedef FieldType_ FieldType;
		
		Phonon(size_t vol,const std::string& mcstarttype) : size(vol),phonon(vol, PsimagLite::Vector<FieldType>(3)),isFrozen(false)
		{
			if (mcstarttype=="none") return;
			IoSimpleIn ioin(mcstarttype);
			(*this)<=ioin;
			if (phonon.size()==0) throw std::runtime_error("Problem in phonon\n");
		}
		
		size_t size;
                //PsimagLite::Vector<FieldType> dummy_;
		std::vector<PsimagLite::Vector<FieldType> > phonon;
		bool isFrozen;
		
	}; // Spin
	
	template<typename FieldType>
	std::ostream& operator<<(std::ostream& os,Phonon<FieldType>& dynVars)
	{
		os<<"Phonon\n";
		os<<dynVars.phonon;
		os<<"IsFrozenPhonon "<<dynVars.isFrozen<<"\n";
		return os;
	}
	
	//! Operator to read Dynvars from file
	template<typename FieldType>
	Phonon<FieldType>&
	operator <= (Phonon<FieldType>& dynVars,  typename Dmrg::IoSimple::In& ioin) 
	{
		ioin.read(dynVars.phonon,"Phonon");
		ioin.readline(dynVars.isFrozen,"IsFrozenPhonon");
		
		return dynVars;
	}

	
} // namespace Spf

/*@}*/
#endif
