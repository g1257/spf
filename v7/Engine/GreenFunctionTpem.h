
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
		typedef std::vector<RealType> VectorType;
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
		  energyCoeffs_(algorithm_.cutoff()),
		  numberCoeffs_(algorithm_.cutoff()),
		  energyFunctor_(algorithm_.tpemParameters()),
		  numberFunctor_(algorithm_.tpemParameters())
		{
			computeCoeffs(energyCoeffs_,energyFunctor_);
			computeCoeffs(numberCoeffs_,numberFunctor_);
		}
		
		void measure()
		{
			algorithm_.prepare();
// 			size_t n = hilbertSize_;
// 			for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++)
// 				data_(i,j) = greenFunction(i,j);
		}

		AlgorithmType& algorithm() { return algorithm_; }

		ModelType& model() { return algorithm_.model(); } // should be const

		size_t hilbertSize() const { return hilbertSize_; }

		const ComplexType& operator()(size_t lambda1,size_t lambda2) const
		{
			throw std::runtime_error("gf is unimplemented for TPEM\n");
// 			return data_(lambda1,lambda2);
		}

		RealType calcNumber() const
		{
// 			tmpValues(aux.varTpem_a,aux.varTpem_b,aux.varMu,beta,0);
// 			if (ether.isSet("adjusttpembounds"))
// 				tpem_calculate_coeffs (numberCoeffs_,numberFunctor_,tpemOptions_);
			
			return tpem_.expand(algorithm_.moment(), numberCoeffs_);
		}

		RealType calcElectronicEnergy() const
		{
// 			tmpValues(aux.varTpem_a,aux.varTpem_b,aux.varMu,beta,0);
// 			if (ether.isSet("adjusttpembounds"))
// 				tpem_calculate_coeffs (energyCoeffs_,energyFunctor,tpemOptions_);
			
			return tpem_.expand(algorithm_.moment(), energyCoeffs_);
		}

		ComplexType matrix(size_t lambda1,size_t lambda2) const
		{
			return ComplexType(0,0);
		}

		FieldType e(size_t i) const
		{
			return 0;
		}

		void localCharge(std::vector<RealType>& lc)
		{
			throw std::runtime_error("local charge is unimplemented for TPEM\n");
// 			//checkUs();
// 			//checkLevels();
// 			for (size_t i=0;i<hilbertSize_;i++) {
// 				for (size_t lambda=0;lambda<hilbertSize_;lambda++) {
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

		void electronSpin(std::vector<ComplexType>& es,
				size_t orbitals,size_t n) const
		{
			throw std::runtime_error("electron spin is unimplemented for TPEM\n");
// 			enum {SPIN_UP,SPIN_DOWN};
// 
// 			for (size_t i=0;i<n;i++) {
// 				for (size_t orb=0;orb<orbitals;orb++) {
// 					size_t x = i+(orb+SPIN_UP*orbitals)*n;
// 					size_t y = i+(orb+SPIN_DOWN*orbitals)*n;
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


		template<typename EngineParametersType2,typename ModelType2,typename RngType2>
		friend std::ostream& operator<<(std::ostream& os,
		  const GreenFunctionTpem<EngineParametersType2,ModelType2,RngType2>& gf);

	private:
		
		template<typename ObservableFunctorType>
		void computeCoeffs(std::vector<RealType>& coeffs,
		                   ObservableFunctorType& obsFunc)
		{
// 			tmpValues(aux.varTpem_a,aux.varTpem_b,aux.varMu,beta,0);
			tpem_.calcCoeffs(coeffs, obsFunc);
		}
		
// 		ComplexType greenFunction(size_t lambda1,size_t lambda2) const
// 		{
// 			ComplexType sum = 0;
// 			RealType beta = engineParams_.beta;
// 			RealType mu = engineParams_.mu;
// 
// 			for (size_t lambda=0;lambda<hilbertSize_;lambda++)
// 				sum += std::conj(algorithm_.matrix(lambda1,lambda)) *
// 					algorithm_.matrix(lambda2,lambda) *
// 				PsimagLite::fermi(-beta*(algorithm_.e(lambda)-mu));
// 			return sum;
// 		}


		const EngineParametersType& engineParams_;		
		AlgorithmType algorithm_;
		size_t hilbertSize_;
		TpemType& tpem_; // we don't own it, FIXME: must be const
		PsimagLite::Matrix<ComplexType> data_;
		VectorType energyCoeffs_;
		VectorType numberCoeffs_;
		EnergyFunctorType energyFunctor_;
		NumberFunctorType numberFunctor_;
		
	}; // class GreenFunctionTpem
	
	template<typename EngineParametersType,typename ModelType,typename RngType>
	std::ostream& operator<<(std::ostream& os,const GreenFunctionTpem<EngineParametersType,ModelType,RngType>& gf)
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
