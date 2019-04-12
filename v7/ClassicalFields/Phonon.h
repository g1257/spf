
/** \ingroup SPF */
/*@{*/

/*! \file Phonon.h
 *
 *  DynVars: classical phonons
 *
 */
#ifndef PHONON_H
#define PHONON_H
#include "Io/IoSimple.h"
#include "Vector.h"
#include "ConstantVector.h"

namespace Spf {

	template<typename FieldType_>
	struct Phonon { // Do not add functions here, this is a struct!!

		typedef typename PsimagLite::IoSimple::In IoSimpleIn;
		typedef FieldType_ FieldType;
		typedef typename PsimagLite::Vector<FieldType>::Type OnePhononType;

		Phonon(SizeType vol,const PsimagLite::String& mcstarttype)
			: size(vol),
		      phonon(vol, typename PsimagLite::Vector<FieldType>::Type(3)),isFrozen(false)
		{
			if (mcstarttype=="none") return;
			IoSimpleIn ioin(mcstarttype);

			ioin.read(phonon,"Phonon");
			ioin.readline(isFrozen,"IsFrozenPhonon");

			if (phonon.size()==0) throw PsimagLite::RuntimeError("Problem in phonon\n");
		}
		
		SizeType size;
                //PsimagLite::Vector<FieldType> dummy_;
		typename PsimagLite::Vector<OnePhononType>::Type phonon;
		ConstantVector modulus;
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
