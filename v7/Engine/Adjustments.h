
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
			
			for (size_t iter=0;iter<maxIter_;iter++) {
				FieldType denom=nOfElectPrime(mu,engineParams_.beta,eigs);
				FieldType tmp = nOfElectrons(mu,engineParams_.beta,eigs);
				if (fabs(denom)<1e-3) break;
				mu = mu -(tmp-n0)/denom;

				if (fabs(nOfElectrons(mu,engineParams_.beta,eigs)-n0)<tolerance_) {
					converged=true;
					break;
				}
			}
			
			if (!converged) {
				std::cerr<<"achieved: "<<nOfElectrons(mu,engineParams_.beta,eigs)<<" "<<mu<<"\n";
				throw std::runtime_error("Adjustments: adjChemPot: Failed to converged.\n");
			}
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

