
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
		typedef typename AlgorithmType::ComplexOrRealType ComplexOrRealType;
		typedef std::complex<RealType> ComplexType;
		typedef RngType RandomNumberGeneratorType;
		typedef typename EngineParametersType::IoInType IoInType;

		GreenFunctionDiag(const EngineParametersType& engineParams,
		                  ModelType& model,
		                  IoInType& io)
		: engineParams_(engineParams),
		  algorithm_(engineParams,model,io),
		  hilbertSize_(model.hilbertSize()),
		  data_(hilbertSize_,hilbertSize_)
		{}

		void measure()
		{
			algorithm_.prepare();
			SizeType n = hilbertSize_;
			for (SizeType i=0;i<n;i++) for (SizeType j=0;j<n;j++)
				data_(i,j) = greenFunction(i,j);
		}

		AlgorithmType& algorithm() { return algorithm_; }

		bool usesDiagonalization() const { return true; }

		ModelType& model() { return algorithm_.model(); } // should be const

		const ComplexOrRealType& operator()(SizeType lambda1,SizeType lambda2) const
		{
			return data_(lambda1,lambda2);
		}

		RealType calcNumber() const
		{
			RealType sum=0;
			for (SizeType i=0;i<hilbertSize_;i++) {
				sum += PsimagLite::fermi(
					(algorithm_.e(i)-engineParams_.mu)*engineParams_.beta);
			}
			return sum;
		}

		RealType calcElectronicEnergy() const
		{
			RealType sum=0;
			for (SizeType i=0;i<hilbertSize_;i++) {
				sum +=algorithm_.e(i)*PsimagLite::fermi(
					(algorithm_.e(i)-engineParams_.mu)*engineParams_.beta);
			}
			return sum;

		}

		void localCharge(typename PsimagLite::Vector<RealType>::Type& lc)
		{
			//checkUs();
			//checkLevels();
			for (SizeType i=0;i<hilbertSize_;i++) {
				for (SizeType lambda=0;lambda<hilbertSize_;lambda++) {
					ComplexOrRealType tmp =conj(algorithm_.matrix(i,lambda))*algorithm_.matrix(i,lambda);
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

		void electronSpin(typename PsimagLite::Vector<ComplexType>::Type& es,
				SizeType orbitals,SizeType n) const
		{
			enum {SPIN_UP,SPIN_DOWN};

			for (SizeType i=0;i<n;i++) {
				for (SizeType orb=0;orb<orbitals;orb++) {
					SizeType x = i+(orb+SPIN_UP*orbitals)*n;
					SizeType y = i+(orb+SPIN_DOWN*orbitals)*n;
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

		const ComplexOrRealType& matrix(SizeType lambda1,SizeType lambda2) const
		{
			return algorithm_.matrix(lambda1,lambda2);
		}

		const RealType& e(SizeType i) const
		{
			return algorithm_.e(i);
		}

		SizeType hilbertSize() const { return algorithm_.hilbertSize(); }

		void printMatrix(SizeType mode) const
		{
			algorithm_.printMatrix(mode);
		}

		template<typename EngineParametersType2,typename ModelType2,typename RngType2>
		friend std::ostream& operator<<(std::ostream& os,
		  const GreenFunctionDiag<EngineParametersType2,ModelType2,RngType2>& gf);

	private:

		ComplexOrRealType greenFunction(SizeType lambda1,SizeType lambda2) const
		{
			ComplexOrRealType sum = 0;
			RealType beta = engineParams_.beta;
			RealType mu = engineParams_.mu;

			for (SizeType lambda=0;lambda<hilbertSize_;lambda++)
				sum += PsimagLite::conj(algorithm_.matrix(lambda1,lambda)) *
					algorithm_.matrix(lambda2,lambda) *
				PsimagLite::fermi(-beta*(algorithm_.e(lambda)-mu));
			return sum;
		}

		void checkUs() const
		{
			for (SizeType i=0;i<hilbertSize_;i++) {
				ComplexOrRealType s = 0;
				for (SizeType lambda=0;lambda<hilbertSize_;lambda++) {
					ComplexOrRealType tmp =conj(algorithm_.matrix(lambda,i))*algorithm_.matrix(lambda,i);
					s += tmp;
				}
				std::cerr<<"i="<<i<<" sum_lambda U*_{lambda,i} U_{lambda,i} = "<<s<<"\n";
			}
		}

		void checkLevels() const
		{
			RealType sum = 0;
			for (SizeType lambda=0;lambda<hilbertSize_;lambda++) {
				if (algorithm_.e(lambda)>=engineParams_.mu) continue;
				sum++;
			}
			std::cerr<<"mu = "<<engineParams_.mu<<" Levels below="<<sum<<"\n";
		}

		const EngineParametersType& engineParams_;

		AlgorithmType algorithm_;
		SizeType hilbertSize_;
		PsimagLite::Matrix<ComplexOrRealType> data_;

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
