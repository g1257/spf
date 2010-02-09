
/** \ingroup SPF */
/*@{*/

/*! \file MonteCarlo.h
 *
 *  Monte Carlo for SPF
 *
 */
#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H
#include "Utils.h"
#include "ProgressIndicator.h"

namespace Spf {
	template<typename EngineParamsType,typename ModelType,typename DynVarsType,typename RandomNumberGeneratorType>
	class MonteCarlo {
		typedef typename DynVarsType::FieldType FieldType;
		
	public:
		MonteCarlo(const EngineParamsType& engineParams,ModelType& model,RandomNumberGeneratorType& rng) 
			: engineParams_(engineParams),model_(model),rng_(rng) { }
		
		size_t operator()(DynVarsType& dynVars, size_t iter)
		{
			size_t acc = 0;
			std::vector<FieldType> eigOld;
			model_.fillAndDiag(eigOld,dynVars);
			size_t n = dynVars.size();
			
			model_.set(dynVars);
			for (size_t i=0;i<n;i++) {
					
				model_.propose(i);
				
				FieldType dsDirect = model_.deltaDirect(i);
				
				//FieldType oldmu=engineParams_.mu;
				std::vector<FieldType> eigNew;
				model_.fillAndDiag(eigNew);
				
				model_.adjustChemPot(eigNew); //changes engineParams_.mu
				FieldType integrationMeasure = model_.integrationMeasure(i);
				
				bool flag=doMetropolis(eigNew,eigOld,dsDirect,integrationMeasure);
					
				if (flag && !dynVars.isFrozen) { // Accepted
					model_.accept(i);
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
		
		const EngineParamsType& engineParams_;
		ModelType& model_;
		RandomNumberGeneratorType& rng_;
		
	}; // MonteCarlo
} // namespace Spf

/*@}*/
#endif // SPF_ENGINE_H

