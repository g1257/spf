
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
		typedef typename GreenFunctionType::FieldType FieldType;
		typedef typename GreenFunctionType::ComplexType ComplexType;
		typedef PsimagLite::Matrix<FieldType> MatrixType;

	public:	

		Conductance(
				const EngineParametersType& engineParams,
				const GreenFunctionType& gf)
		: engineParams_(engineParams),gf_(gf)
		{
		}

		// needs e^2 \pi
		FieldType operator()(const MatrixType& velocitySquared) const
		{
			FieldType sum = 0;
			FieldType beta = engineParams_.beta;
			FieldType mu = engineParams_.mu;
			FieldType eps = 0.001;
			for (size_t a=0;a<gf_.hilbertSize();a++) {
				FieldType ea = PsimagLite::fermi(beta*(gf_.e(a)-mu));
				FieldType eMa = PsimagLite::fermi(-beta*(gf_.e(a)-mu));
				for (size_t b=0;b<gf_.hilbertSize();b++) {
					if (fabs(gf_.e(a)-gf_.e(b))>eps) continue;
					//FieldType eb = PsimagLite::fermi(beta*(gf_.e(b)-mu));
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
