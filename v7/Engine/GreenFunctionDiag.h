
/** \ingroup SPF */
/*@{*/

/*! \file GreenFunctionDiag.h
 *
 *  FIXME
 *
 */
#ifndef GREENFUNCTION_DIAG_H
#define GREENFUNCTION_DIAG_H
#include "Matrix.h" // in PsimagLite
#include "Fermi.h" // in PsimagLite
#include "AlgorithmDiag.h"

namespace Spf {
	template<typename EngineParametersType,typename ModelType_,typename RngType>
	class GreenFunctionDiag {
	public:
		typedef ModelType_ ModelType;
		typedef AlgorithmDiag<EngineParametersType,ModelType_,RngType>
		AlgorithmType;
		typedef typename AlgorithmType::RealType RealType;
		typedef typename AlgorithmType::ComplexType ComplexType;
		typedef RngType RandomNumberGeneratorType;

		GreenFunctionDiag(const EngineParametersType& engineParams,ModelType& model)
		: engineParams_(engineParams),
		  algorithm_(engineParams,model),
		  hilbertSize_(model.hilbertSize()),
		  data_(hilbertSize_,hilbertSize_)
		{}
		
		void measure()
		{
			algorithm_.prepare();
			size_t n = hilbertSize_;
			for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++)
				data_(i,j) = greenFunction(i,j);
		}
		
		AlgorithmType& algorithm() { return algorithm_; }
		
		ModelType& model() { return algorithm_.model(); } // should be const
		
		const ComplexType& operator()(size_t lambda1,size_t lambda2) const
		{
			return data_(lambda1,lambda2);
		}
		
		RealType calcNumber() const
		{
			RealType sum=0;
			for (size_t i=0;i<hilbertSize_;i++) {
				sum += PsimagLite::fermi(
					(algorithm_.e(i)-engineParams_.mu)*engineParams_.beta);
			}
			return sum;
		}
		
		RealType calcElectronicEnergy() const
		{
			RealType sum=0;
			for (size_t i=0;i<hilbertSize_;i++) {
				sum +=algorithm_.e(i)*PsimagLite::fermi(
					(algorithm_.e(i)-engineParams_.mu)*engineParams_.beta);
			}
			return sum;
				
		}
		
		void localCharge(std::vector<RealType>& lc)
		{
			//checkUs();
			//checkLevels();
			for (size_t i=0;i<hilbertSize_;i++) {
				for (size_t lambda=0;lambda<hilbertSize_;lambda++) {
					ComplexType tmp =conj(algorithm_.matrix(i,lambda))*algorithm_.matrix(i,lambda);
					//if (algorithm_.e(lambda)>=engineParams_.mu) continue; // temperature zero
					RealType s = real(tmp)*PsimagLite::fermi(engineParams_.beta*
							(algorithm_.e(lambda)-engineParams_.mu));
					//if (ether.isSet("savelcd")) {
					//	lc[i+alpha*linSize] = s;
					//} else {
						lc[i] += s;
					//}
				}
				//std::cerr<<"Done i="<<i<<" ss="<<ss<<"\n";
			}
		}

		void electronSpin(std::vector<ComplexType>& es,
				size_t orbitals,size_t n) const
		{
			enum {SPIN_UP,SPIN_DOWN};

			for (size_t i=0;i<n;i++) {
				for (size_t orb=0;orb<orbitals;orb++) {
					size_t x = i+(orb+SPIN_UP*orbitals)*n;
					size_t y = i+(orb+SPIN_DOWN*orbitals)*n;
					// x-component
					es[0] -= data_(x,y) + data_(y,x);
					// y-component
					es[1] += data_(x,y) - data_(y,x);
					// z-component
					es[2] += (1.0 - data_(x,x));
					es[2] -= (1.0 - data_(y,y));
				}
			}
			// y component needs multiplication by sqrt(-1):
			es[1] *= ComplexType(0,1);
		}

		const ComplexType& matrix(size_t lambda1,size_t lambda2) const
		{
			return algorithm_.matrix(lambda1,lambda2);
		}

		const RealType& e(size_t i) const
		{
			return algorithm_.e(i);
		}

		size_t hilbertSize() const { return algorithm_.hilbertSize(); }
		
		void printMatrix(size_t mode) const
		{
			algorithm_.printMatrix(mode);
		}

		template<typename EngineParametersType2,typename ModelType2,typename RngType2>
		friend std::ostream& operator<<(std::ostream& os,
		  const GreenFunctionDiag<EngineParametersType2,ModelType2,RngType2>& gf);

	private:
		
		ComplexType greenFunction(size_t lambda1,size_t lambda2) const
		{
			ComplexType sum = 0;
			RealType beta = engineParams_.beta;
			RealType mu = engineParams_.mu;

			for (size_t lambda=0;lambda<hilbertSize_;lambda++)
				sum += std::conj(algorithm_.matrix(lambda1,lambda)) *
					algorithm_.matrix(lambda2,lambda) *
				PsimagLite::fermi(-beta*(algorithm_.e(lambda)-mu));
			return sum;
		}

		void checkUs() const
		{
			for (size_t i=0;i<hilbertSize_;i++) {
				ComplexType s = 0;
				for (size_t lambda=0;lambda<hilbertSize_;lambda++) {
					ComplexType tmp =conj(algorithm_.matrix(lambda,i))*algorithm_.matrix(lambda,i);
					s += tmp;
				}
				std::cerr<<"i="<<i<<" sum_lambda U*_{lambda,i} U_{lambda,i} = "<<s<<"\n";
			}
		}

		void checkLevels() const
		{
			RealType sum = 0;
			for (size_t lambda=0;lambda<hilbertSize_;lambda++) {
				if (algorithm_.e(lambda)>=engineParams_.mu) continue;
				sum++;
			}
			std::cerr<<"mu = "<<engineParams_.mu<<" Levels below="<<sum<<"\n"; 
		}

		const EngineParametersType& engineParams_;
		
		AlgorithmType algorithm_;
		size_t hilbertSize_;
		PsimagLite::Matrix<ComplexType> data_;
		
	}; // GreenFunctionDiag
	
	template<typename EngineParametersType,typename ModelType,typename RngType>
	std::ostream& operator<<(std::ostream& os,const GreenFunctionDiag<EngineParametersType,ModelType,RngType>& gf)
	{
		os<<"#GF\n";
		os<<gf.data_;
		os<<gf.algorithm_;
		return os;
	}
} // namespace Spf

/*@}*/
#endif
