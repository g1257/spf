
/** \ingroup SPF */
/*@{*/

/*! \file MonteCarloHelper.h
 *
 *  Monte Carlo Helper for SPF
 *
 */
#ifndef MONTE_CARLO_HELPER_H
#define MONTE_CARLO_HELPER_H
#include "Utils.h"
#include "ProgressIndicator.h"

namespace Spf {
	template<typename MonteCarloHelper,typename DynVarsType>
	class MonteCarlo {
	public:
		size_t operator()(DynVarsType& dynVars, MonteCarloHelper& mcHelper, size_t iter)
		{
			size_t acc = 0;
			std::vector<FieldType> eigOld;
			mcHelper.fillAndDiag(eigOld,dynVars);
			size_t n = dynVars.size();
			
			operations.set(dynVars);
			for (size_t i=0;i<n;i++) {
					
				mcHelper.propose(i);
				
				FieldType dsDirect = mcHelper.deltaDirect(i);
				
				//FieldType oldmu=engineParams_.mu;
				std::vector<FieldType> eigNew;
				mcHelper.fillAndDiag(eigNew,classicalSpinOperations_.dynVars2());
				
				adjustChemPot(eigNew); //changes engineParams_.mu
				FieldType integrationMeasure = mcHelper.integrationMeasure(i);
				
				bool flag=doMetropolis(eigNew,eigOld,dsDirect,integrationMeasure);
					
				if (flag && !dynVars.isFrozen) { // Accepted
					operations.accept(i);
					eigOld = eigNew;
					acc++;
				} else { // not accepted
					//engineParams_.mu=oldmu;
				}
			} // lattice sweep
			return acc;
		}
	private:
		bool doMetropolis(const std::vector<FieldType>& eNew,const std::vector<FieldType>& eOld,
			FieldType dsDirect,FieldType integrationMeasure)
		{
			FieldType mu=engineParams_.mu;
			FieldType beta = engineParams_.beta;
			FieldType X =1.0;
			
			for (size_t i=0;i<eNew.size();i++) {
				FieldType temp = 0;
				if (eNew[i]>mu)
					temp = (double)(1.0+exp(-beta*(eNew[i]-mu)))/(1.0+exp(-beta*(eOld[i]-mu)));
				else
				temp =(double)(1.0+exp(beta*(eNew[i]-mu)))/
							(exp(beta*(eNew[i]-mu))+exp(-beta*(eOld[i]-eNew[i])));
			
				X *= temp;
			}
			
			//if (ether.isSet("sineupdate")) X *= integrationMeasure;
			X *=  exp(-beta*dsDirect);
			X = X/(1.0+X);
			
			FieldType r=rng_();
			
			if (X>r) return true;
			else return false;
		}	
	}; // MonteCarlo
} // namespace Spf

/*@}*/
#endif // SPF_ENGINE_H

