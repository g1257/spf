
/** \ingroup SPF */
/*@{*/

/*! \file RootFindingNewton.h
 *
 *  RootFindingNewton such as chemical potential
 *
 */
#ifndef ROOT_FIND_NEWTON_H
#define ROOT_FIND_NEWTON_H
#include "ProgressIndicator.h"
#include "Matrix.h" // in PsimagLite
#include "Fermi.h"// in PsimagLite

/* ****************** ALIAGA'S ALGORITHM ****************************/
// ADJUST CHEMICAL POTENTIAL FOR DIAGONALIZATION ****************/

namespace Spf {
	template<typename FunctionType>
	class RootFindingNewton {
		typedef typename FunctionType::FieldType FieldType;
		
	public:
		RootFindingNewton(const FunctionType& function,
		                  SizeType maxIter=100,
		                  FieldType tolerance=1.0e-3)
		: function_(function),maxIter_(maxIter),tolerance_(tolerance)
		{
		}
		
		void operator()(FieldType& mu) const
		{
			bool converged=false;
			//std::cerr<<"Starting with "<<mu<<" n0="<<n0<<"\n";

			for (SizeType iter=0;iter<maxIter_;iter++) {
				FieldType denom=function_.derivative(mu);
				FieldType tmp = function_(mu);
				//std::cerr<<"with mu="<<mu<<" denom(prime) = "<<denom<<" electrons=" <<tmp<<"\n";

				//std::cerr<<"denom="<<denom<<"\n";
				if (fabs(denom)<1e-5) break;
				mu = mu -(tmp)/denom;
				//std::cerr<<"iter="<<iter<<" mu ="<<mu<<"\n";
				if (fabs(function_(mu))<tolerance_) {
					converged=true;
					break;
				}
			}
			FieldType nn = function_(mu);
			if (fabs(nn)<tolerance_) {
				converged=true;
			}

			if (!converged) {
				std::cerr<<"achieved: "<<nn<<" "<<mu<<"\n";
				throw std::runtime_error(
				      "RootFindingNewton: adjChemPot: Failed to converged.\n");
			}
			//std::cerr<<"Converged: "<<mu<<" "<<nn<<"\n";
		}
		
	private:
		
		const FunctionType& function_;
		SizeType maxIter_;
		FieldType tolerance_;
	}; // RootFindingNewton
} // namespace Spf

/*@}*/
#endif// ROOT_FIND_NEWTON_H

