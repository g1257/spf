
/** \ingroup SPF */
/*@{*/

/*! \file Conductance.h
 *
 *  DOCUMENTATION NEEDED FIXME
 *
 */
#ifndef CONDUCTANCE_H
#define CONDUCTANCE_H
#include "Fermi.h" // in PsimagLite

namespace Spf {
	template<typename EngineParametersType,typename GreenFunctionType>
	class Conductance {
		typedef typename GreenFunctionType::RealType RealType;
		typedef typename GreenFunctionType::ComplexType ComplexType;
		typedef PsimagLite::Matrix<RealType> MatrixType;

	public:	

		Conductance(
				const EngineParametersType& engineParams,
				const GreenFunctionType& gf)
		: engineParams_(engineParams),gf_(gf)
		{
		}

		// needs e^2 \pi
		RealType operator()(const MatrixType& velocitySquared) const
		{
			RealType sum = 0;
			RealType beta = engineParams_.beta;
			RealType mu = engineParams_.mu;
			RealType eps = 0.001;
			for (size_t a=0;a<gf_.hilbertSize();a++) {
				RealType ea = PsimagLite::fermi(beta*(gf_.e(a)-mu));
				RealType eMa = PsimagLite::fermi(-beta*(gf_.e(a)-mu));
				for (size_t b=0;b<gf_.hilbertSize();b++) {
					if (fabs(gf_.e(a)-gf_.e(b))>eps) continue;
					//RealType eb = PsimagLite::fermi(beta*(gf_.e(b)-mu));
// 					if (fabs(ea-eb)>eps) continue;
					sum += velocitySquared(a,b)*beta*eMa*ea;
				}
			}
			return sum;
		}

		const EngineParametersType& engineParams_;	
		const GreenFunctionType& gf_;
		
	}; // Conductance
	
//	template<typename EngineParametersType,typename AlgorithmType>
//	std::ostream& operator<<(std::ostream& os,Conductance<EngineParametersType,AlgorithmType>& algorithm)
//	{
//		throw std::runtime_error("unimplemented operator<< for Conductance\n");
//		return os;
//	}
} // namespace Spf

/*@}*/
#endif //CONDUCTANCE_H
