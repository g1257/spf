
/** \ingroup SPF */
/*@{*/

/*! \file Adjustments.h
 *
 *  Adjustments such as chemical potential
 *
 */
#ifndef ADJUSTMENTS_H
#define ADJUSTMENTS_H
#include "Utils.h"
#include "ProgressIndicator.h"

/* ****************** ALIAGA'S ALGORITHM ****************************/
// ADJUST CHEMICAL POTENTIAL FOR DIAGONALIZATION ****************/

namespace Spf {
	template<typename EngineParamsType>
	class Adjustments {
		typedef typename EngineParamsType::FieldType FieldType;
		
	public:
		Adjustments(const EngineParamsType& engineParams,size_t maxIter=100,FieldType tolerance=1.0e-3) :
			engineParams_(engineParams),maxIter_(maxIter),tolerance_(tolerance)
		{
		}
		
		FieldType adjChemPot(const std::vector<FieldType>& eigs) const
		{
			FieldType mu=engineParams_.mu;
			bool converged=false;
			size_t n0=engineParams_.carriers;
			//std::cerr<<"Starting with "<<mu<<" n0="<<n0<<"\n";

			for (size_t iter=0;iter<maxIter_;iter++) {
				FieldType denom=nOfElectPrime(mu,engineParams_.beta,eigs);
				FieldType tmp = nOfElectrons(mu,engineParams_.beta,eigs);
				//std::cerr<<"with mu="<<mu<<" denom(prime) = "<<denom<<" electrons=" <<tmp<<"\n";

				//std::cerr<<"denom="<<denom<<"\n";
				if (fabs(denom)<1e-5) break;
				mu = mu -(tmp-n0)/denom;
				//std::cerr<<"iter="<<iter<<" mu ="<<mu<<"\n";
				if (fabs(nOfElectrons(mu,engineParams_.beta,eigs)-n0)<tolerance_) {
					converged=true;
					break;
				}
			}
			FieldType nn = nOfElectrons(mu,engineParams_.beta,eigs);
			if (fabs(nn-n0)<tolerance_) {
				converged=true;
			}

			if (!converged) {
				std::cerr<<"achieved: "<<nn<<" "<<mu<<"\n";
				throw std::runtime_error("Adjustments: adjChemPot: Failed to converged.\n");
			}
			//std::cerr<<"Converged: "<<mu<<" "<<nn<<"\n";
			return mu;
		}
		
		void print(std::ostream& os) const
		{
			os<<"Adjustments: mu="<<engineParams_.mu<<"\n";
		}
		
	private:

		FieldType nOfElectrons(FieldType mu,FieldType beta,const std::vector<FieldType>& eig) const
		{
			FieldType sum = 0;
			for (size_t i=0;i<eig.size();i++) sum += utils::fermi((eig[i]-mu)*beta);
			//for (size_t i=0;i<eig.size()-1;i++) {
			       //if (eig[i]>eig[i+1]) throw std::runtime_error("problem sorting\n");
			       //std::cerr<<i<<" "<<eig[i]<<"\n";
			//}	       
			return sum;
		}

		// Derivative of n(mu) with respect to mu
		FieldType nOfElectPrime(FieldType mu,FieldType beta,const std::vector<FieldType>& eig) const
		{
			FieldType sum=0;
			for (size_t i=0;i<eig.size();i++) sum -= utils::fermiPrime((eig[i]-mu)*beta)*beta;
			return sum;
		}
		
		const EngineParamsType& engineParams_;
		size_t maxIter_;
		FieldType tolerance_;
	}; // Adjustments
} // namespace Spf

/*@}*/
#endif// ADJUSTMENTS_H

