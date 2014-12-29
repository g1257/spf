
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
#include "Random48.h"
#include "TypeToString.h"

namespace Spf {
	template<typename FieldType_>
	struct Spin { // Do not add functions here, this is a struct!!
		typedef typename PsimagLite::IoSimple::In IoSimpleIn;
		typedef FieldType_ FieldType;

		template<typename SomeParamsType>
		Spin(SizeType vol,const SomeParamsType& params)
		    : size(vol),
		      theta(vol,0),
		      phi(vol,0),
		      modulus(vol,1),
		      isFrozen(params.options.find("frozenspins") != PsimagLite::String::npos)
		{
			if (params.dynvarsfile=="none") return;
			if (params.dynvarsfile==":random") {
				PsimagLite::Random48<FieldType> random48(params.randomSeed);
				for (SizeType i=0;i<theta.size();i++) {
					theta[i] = random48.random()*M_PI;
					phi[i] = random48.random()*2.*M_PI;
				}
				return;
			}

			if (params.dynvarsfile==":fm") {
				for (SizeType i=0;i<theta.size();i++) {
					theta[i] = 0.0;
					phi[i] = 0.0;
				}
				return;
			}

			if (params.dynvarsfile==":pizero" || params.dynvarsfile==":afm") {
				SizeType l = SizeType(sqrt(size));
				if (l*l!=size) {
					PsimagLite::String s(__FILE__);
					s += " : " + ttos(__LINE__);
					s += ": Hi there, I'm the Spin class, I have no way of ";
					s += " knowing what geometry you are using, but ";
					s += " it sure doesn't appear to be a square lattice.\n";
					s += " \"pizero\" start type valid only for square lattice\n";
					throw PsimagLite::RuntimeError(s.c_str());
				}

				if (params.dynvarsfile==":afm") {
					for (SizeType i=0;i<theta.size();i++) {
						div_t q = div(i,l);
						theta[i] = ((q.quot+q.rem)&1) ? 0.0 : M_PI;
						phi[i] = 0.0;
					}

				} else {
					for (SizeType i=0;i<theta.size();i++) {
						theta[i] = M_PI;
						phi[i] = (i % l) % 2 ? 0 : M_PI;
					}
				}

				return;
			}

			IoSimpleIn ioin(params.dynvarsfile);
			(*this)<=ioin;
			if (theta.size()==0 || phi.size()==0) throw PsimagLite::RuntimeError("PRoblem\n");
		}

		SizeType size;
		typename PsimagLite::Vector<FieldType>::Type theta;
		typename PsimagLite::Vector<FieldType>::Type phi;
		typename PsimagLite::Vector<SizeType>::Type modulus;
		bool isFrozen;
	}; // Spin

	template<typename FieldType>
	std::ostream& operator<<(std::ostream& os,const Spin<FieldType>& dynVars)
	{
		os<<"Theta\n";
		os<<dynVars.theta;
		os<<"Phi\n";
		os<<dynVars.phi;
		os<<"modulus\n";
		for (SizeType i=0;i<dynVars.modulus.size();i++)
			if (dynVars.modulus[i]!=0) os<<i<<" ";
		os<<"\n";
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
