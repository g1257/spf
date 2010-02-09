
/** \ingroup SPF */
/*@{*/

/*! \file AlgorithmDiag.h
 *
 *  Diagonalization method for SPF
 *
 */
#ifndef ALGORITHM_DIAG_H
#define ALGORITHM_DIAG_H
#include "Utils.h"
#include "ProgressIndicator.h"

namespace Spf {
	template<typename EngineParametersType,typename ModelType,typename RandomNumberGeneratorType>
	class AlgorithmDiag {
		
		typedef typename EngineParametersType::FieldType FieldType;
	public:
		
		AlgorithmDiag(const EngineParametersType& engineParams,ModelType& model)
			: engineParams_(engineParams),model_(model),rng_(),
					eigNew_(model.hilbertSize()),eigOld_(model.hilbertSize())
		{
		}
		
		
		template<typename DynVarsType>
		void init(DynVarsType& dynVars)
		{
			model_.fillAndDiag(eigOld_,dynVars);
			model_.set(dynVars);
		}	
		
		bool isAccepted(size_t i)
		{
			FieldType dsDirect = model_.deltaDirect(i);
				
			//FieldType oldmu=engineParams_.mu;
			
			model_.fillAndDiag(eigNew_);
				
			model_.adjustChemPot(eigNew_); //changes engineParams_.mu
			FieldType integrationMeasure = model_.integrationMeasure(i);
				
			return doMetropolis(dsDirect,integrationMeasure);
		}
		
		void accept(size_t i)
		{
			model_.accept(i);
			eigOld_ = eigNew_;
		}
		
	private:
		bool doMetropolis(FieldType dsDirect,FieldType integrationMeasure)
		{
			FieldType mu=engineParams_.mu;
			FieldType beta = engineParams_.beta;
			FieldType X =1.0;
			
			for (size_t i=0;i<eigNew_.size();i++) {
				FieldType temp = 0;
				if (eigNew_[i]>mu)
					temp = (double)(1.0+exp(-beta*(eigNew_[i]-mu)))/(1.0+exp(-beta*(eigOld_[i]-mu)));
				else
				temp =(double)(1.0+exp(beta*(eigNew_[i]-mu)))/
							(exp(beta*(eigNew_[i]-mu))+exp(-beta*(eigOld_[i]-eigNew_[i])));
			
				X *= temp;
			}
			
			//if (ether.isSet("sineupdate")) X *= integrationMeasure;
			X *=  exp(-beta*dsDirect);
			X = X/(1.0+X);
			
			FieldType r=rng_();
			
			if (X>r) return true;
			else return false;
		}	
		
		const EngineParametersType& engineParams_;
		ModelType& model_;
		RandomNumberGeneratorType rng_;
		std::vector<FieldType> eigNew_,eigOld_;
	}; // AlgorithmDiag
} // namespace Spf

/*@}*/
#endif // ALGORITHM_DIAG_H
