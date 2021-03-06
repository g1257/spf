
/** \ingroup SPF */
/*@{*/

/*! \file GreenFunctionTpem.h
 *
 *  FIXME
 *
 */
#ifndef GREENFUNCTION_TPEM_H
#define GREENFUNCTION_TPEM_H
#include "Matrix.h" // in PsimagLite
#include "Fermi.h" // in PsimagLite
#include "AlgorithmTpem.h"
#include <cmath>

namespace Spf {
	template<typename EngineParametersType,typename ModelType_,typename RngType>
	class GreenFunctionTpem {
	public:
		typedef ModelType_ ModelType;
		typedef AlgorithmTpem<EngineParametersType,ModelType_,RngType>
		AlgorithmType;
		typedef typename AlgorithmType::RealType RealType;
		typedef typename AlgorithmType::ComplexType ComplexType;
		typedef typename AlgorithmType::TpemType TpemType;
		typedef RngType RandomNumberGeneratorType;
		typedef typename PsimagLite::Vector<RealType>::Type VectorType;
		typedef typename TpemType::EnergyFunctorType EnergyFunctorType;
		typedef typename TpemType::NumberFunctorType NumberFunctorType;
		typedef typename EngineParametersType::IoInType IoInType;

		GreenFunctionTpem(const EngineParametersType& engineParams,
		                  ModelType& model,
		                  IoInType& io)
		: engineParams_(engineParams),
		  algorithm_(engineParams,model,io),
		  hilbertSize_(model.hilbertSize()),
		  tpem_(algorithm_.tpem()),
		  data_(hilbertSize_,hilbertSize_),
		  energyFunctor_(algorithm_.tpemParameters()),
		  numberFunctor_(algorithm_.tpemParameters())
		{}

		void measure()
		{
			algorithm_.prepare();
			data_(0,0) = -1; // unimplemented
// 			SizeType n = hilbertSize_;
// 			for (SizeType i=0;i<n;i++) for (SizeType j=0;j<n;j++)
// 				data_(i,j) = greenFunction(i,j);
		}

		AlgorithmType& algorithm() { return algorithm_; }

		ModelType& model() { return algorithm_.model(); } // should be const

		SizeType hilbertSize() const { return hilbertSize_; }

		bool usesDiagonalization() const { return false; }

		const ComplexType& operator()(SizeType,SizeType) const
		{
			//! gf is unimplemented for TPEM
			return data_(0,0);
// 			return data_(lambda1,lambda2);
		}

		RealType calcNumber() const
		{
// 			tmpValues(aux.varTpem_a,aux.varTpem_b,aux.varMu,beta,0);
// 			if (ether.isSet("adjusttpembounds"))
// 				tpem_calculate_coeffs (numberCoeffs_,numberFunctor_,tpemOptions_);
			VectorType numberCoeffs(algorithm_.tpemParameters().cutoff);
			tpem_.calcCoeffs(numberCoeffs,numberFunctor_);
			RealType tmp =  tpem_.expand(algorithm_.moment(), numberCoeffs);
			return tmp;
		}

		RealType calcElectronicEnergy() const
		{
// 			tmpValues(aux.varTpem_a,aux.varTpem_b,aux.varMu,beta,0);
// 			if (ether.isSet("adjusttpembounds"))
// 				tpem_calculate_coeffs (energyCoeffs_,energyFunctor,tpemOptions_);
// 			SizeType min=5;
// 			for (SizeType i=min;i<energyCoeffs_.size();i++)
// 				std::cerr<<"cutoff="<<i<<" value="<<tpem_.expand(algorithm_.moment(),energyCoeffs_,i)<<"\n";
			VectorType energyCoeffs(algorithm_.tpemParameters().cutoff);
			tpem_.calcCoeffs(energyCoeffs,energyFunctor_);
			return tpem_.expand(algorithm_.moment(), energyCoeffs);
		}

		ComplexType matrix(SizeType,SizeType) const
		{
			return ComplexType(0,0);
		}

		RealType e(SizeType) const
		{
			return 0;
		}

		void localCharge(typename PsimagLite::Vector<RealType>::Type& lc)
		{
// 			throw PsimagLite::RuntimeError("local charge is unimplemented for TPEM\n");
// 			//checkUs();
// 			//checkLevels();
			for (SizeType i=0;i<hilbertSize_;i++) lc[i] = -1;
// 				for (SizeType lambda=0;lambda<hilbertSize_;lambda++) {
// 					ComplexType tmp =conj(algorithm_.matrix(i,lambda))*algorithm_.matrix(i,lambda);
// 					//if (algorithm_.e(lambda)>=engineParams_.mu) continue; // temperature zero
// 					RealType s = real(tmp)*PsimagLite::fermi(engineParams_.beta*
// 							(algorithm_.e(lambda)-engineParams_.mu));
// 					//if (ether.isSet("savelcd")) {
// 					//	lc[i+alpha*linSize] = s;
// 					//} else {
// 						lc[i] += s;
// 					//}
// 				}
// 				//std::cerr<<"Done i="<<i<<" ss="<<ss<<"\n";
// 			}
		}

		void electronSpin(typename PsimagLite::Vector<ComplexType>::Type& es,
				SizeType,SizeType) const
		{
			for (SizeType i=0;i<es.size();i++)  es[i]=0;
			//throw PsimagLite::RuntimeError("electron spin is unimplemented for TPEM\n");
// 			enum {SPIN_UP,SPIN_DOWN};
//
// 			for (SizeType i=0;i<n;i++) {
// 				for (SizeType orb=0;orb<orbitals;orb++) {
// 					SizeType x = i+(orb+SPIN_UP*orbitals)*n;
// 					SizeType y = i+(orb+SPIN_DOWN*orbitals)*n;
// 					// x-component
// 					es[0] -= data_(x,y) + data_(y,x);
// 					// y-component
// 					es[1] += data_(x,y) - data_(y,x);
// 					// z-component
// 					es[2] += (1.0 - data_(x,x));
// 					es[2] -= (1.0 - data_(y,y));
// 				}
// 			}
// 			// y component needs multiplication by sqrt(-1):
// 			es[1] *= ComplexType(0,1);
		}

		PsimagLite::String error() const
		{
			return algorithm_.error();
		}

		template<typename EngineParametersType2,typename ModelType2,
		         typename RngType2>
		friend std::ostream& operator<<(std::ostream& os,
		  const GreenFunctionTpem<EngineParametersType2,ModelType2,
		                          RngType2>& gf);

	private:

		const EngineParametersType& engineParams_;
		AlgorithmType algorithm_;
		SizeType hilbertSize_;
		TpemType& tpem_; // we don't own it, FIXME: must be const
		PsimagLite::Matrix<ComplexType> data_;
		EnergyFunctorType energyFunctor_;
		NumberFunctorType numberFunctor_;

	}; // class GreenFunctionTpem

	template<typename EngineParametersType,typename ModelType,typename RngType>
	std::ostream& operator<<(std::ostream& os,const GreenFunctionTpem<
	              EngineParametersType,ModelType,RngType>& gf)
	{
		os<<"#GF\n";
		os<<"#GF is unimplemented for TPEM\n";
// 		os<<gf.data_;
		os<<gf.algorithm_;
		return os;
	}
} // namespace Spf

/*@}*/
#endif
