
/** \ingroup SPF */
/*@{*/

/*! \file Phonon.h
 *
 *  DynVars: classical phonons
 *
 */
#ifndef PHONON_H
#define PHONON_H
#include "IoSimple.h"
#include "Vector.h"

namespace Spf {
	template<typename FieldType_>
	struct Phonon { // Do not add functions here, this is a struct!!
		typedef typename PsimagLite::IoSimple::In IoSimpleIn;
		typedef FieldType_ FieldType;
		
		Phonon(size_t vol,const PsimagLite::String& mcstarttype)
			: size(vol),phonon(vol, typename PsimagLite::Vector<FieldType>::Type(3)),isFrozen(false)
		{
			if (mcstarttype=="none") return;
			IoSimpleIn ioin(mcstarttype);

			ioin.read(phonon,"Phonon");
			ioin.readline(isFrozen,"IsFrozenPhonon");

			if (phonon.size()==0) throw PsimagLite::RuntimeError("Problem in phonon\n");
		}
		
		size_t size;
                //PsimagLite::Vector<FieldType> dummy_;
		typename PsimagLite::Vector<typename PsimagLite::Vector<FieldType>::Type>::Type phonon;
		bool isFrozen;
		
	}; // Spin
	
	template<typename FieldType>
	std::ostream& operator<<(std::ostream& os,const Phonon<FieldType>& dynVars)
	{
		os<<"Phonon\n";
		os<<dynVars.phonon;
		os<<"IsFrozenPhonon="<<dynVars.isFrozen<<"\n";
		return os;
	}
	
} // namespace Spf

/*@}*/
#endif
