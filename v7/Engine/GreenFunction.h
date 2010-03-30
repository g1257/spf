
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
			engineParams_(engineParams),algorithm_(algorithm),hilbertSize_(hilbertSize)
		{
		}
		
		ComplexType operator()(size_t lambda1,size_t lambda2)
		{
			return algorithm_.greenFunction(lambda1,lambda2);
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
			for (size_t i=0;i<hilbertSize_;i++) {
				for (size_t lambda=0;lambda<hilbertSize_;lambda++) {
					ComplexType tmp =conj(algorithm_.matrix(lambda,i))*algorithm_.matrix(lambda,i);
					FieldType s = real(tmp)*utils::fermi(engineParams_.beta*
							(algorithm_.e(lambda)-engineParams_.mu));
					//if (ether.isSet("savelcd")) {
					//	lc[i+alpha*linSize] = s;
					//} else {
						lc[i] += s;
					//}
				}
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
		const EngineParametersType& engineParams_;	
		AlgorithmType& algorithm_;
		size_t hilbertSize_;
		
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
