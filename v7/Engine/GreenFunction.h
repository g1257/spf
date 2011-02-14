
/** \ingroup SPF */
/*@{*/

/*! \file GreenFunction.h
 *
 *  FIXME
 *
 */
#ifndef GREENFUNCTION_H
#define GREENFUNCTION_H
#include "Utils.h"

namespace Spf {
	template<typename EngineParametersType,typename AlgorithmType>
	class GreenFunction {
		typedef typename AlgorithmType::FieldType FieldType;
		typedef typename AlgorithmType::ComplexType ComplexType;
	public:	
		
		GreenFunction(const EngineParametersType& engineParams,AlgorithmType& algorithm,size_t hilbertSize) :
			engineParams_(engineParams),algorithm_(algorithm),hilbertSize_(hilbertSize),data_(hilbertSize,hilbertSize)
		{
			algorithm_.prepare();
			size_t n = hilbertSize_;
			for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++)
				data_(i,j) = algorithm_.greenFunction(i,j);
		}
		
		ComplexType operator()(size_t lambda1,size_t lambda2)
		{
			return data_(lambda1,lambda2);
		}
		
		FieldType calcNumber() const
		{
			FieldType sum=0;
			for (size_t i=0;i<hilbertSize_;i++) {
				sum += utils::fermi((algorithm_.e(i)-engineParams_.mu)*engineParams_.beta);
			}
			return sum;
		}
		
		FieldType calcElectronicEnergy() const
		{
			FieldType sum=0;
			for (size_t i=0;i<hilbertSize_;i++) {
				sum +=algorithm_.e(i)*utils::fermi((algorithm_.e(i)-engineParams_.mu)*engineParams_.beta);
			}
			return sum;
				
		}
		
		void localCharge(PsimagLite::Vector<FieldType>& lc)
		{
			//checkUs();
			//checkLevels();
			for (size_t i=0;i<hilbertSize_;i++) {
				for (size_t lambda=0;lambda<hilbertSize_;lambda++) {
					ComplexType tmp =conj(algorithm_.matrix(i,lambda))*algorithm_.matrix(i,lambda);
					//if (algorithm_.e(lambda)>=engineParams_.mu) continue; // temperature zero
					FieldType s = real(tmp)*utils::fermi(engineParams_.beta*
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

		ComplexType matrix(size_t lambda1,size_t lambda2)
		{
			return algorithm_.matrix(lambda1,lambda2);
		}
		
		FieldType e(size_t i)
		{
			return algorithm_.e(i);
		}
		
	private:
		
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
			FieldType sum = 0;
			for (size_t lambda=0;lambda<hilbertSize_;lambda++) {
				if (algorithm_.e(lambda)>=engineParams_.mu) continue;
				sum++;
			}
			std::cerr<<"mu = "<<engineParams_.mu<<" Levels below="<<sum<<"\n"; 
		}

		const EngineParametersType& engineParams_;	
		AlgorithmType& algorithm_;
		size_t hilbertSize_;
		psimag::Matrix<ComplexType> data_;
		
	}; // GreenFunction
	
	template<typename EngineParametersType,typename AlgorithmType>
	std::ostream& operator<<(std::ostream& os,GreenFunction<EngineParametersType,AlgorithmType>& algorithm)
	{
		throw std::runtime_error("unimplemented operator<< for GreenFunction\n");
		return os;
	}
} // namespace Spf

/*@}*/
#endif
