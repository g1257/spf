
/** \ingroup SPF */
/*@{*/

/*! \file MonteCarlo.h
 *
 *  Monte Carlo for SPF
 *
 */
#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H
#include "ProgressIndicator.h"
#include "loki/Typelist.h"

namespace Spf {
	template<typename EngineParamsType,
             typename OperationsList,
             typename AlgorithmType,
             typename RandomNumberGeneratorType,
             int n>
	class MonteCarlo {

		typedef typename Loki::TL::TypeAt<OperationsList,n>::Result OperationsType;
		typedef typename OperationsType::DynVarsType DynVarsType;
		typedef std::pair<size_t,size_t> PairType;
		
	public:
		typedef typename EngineParamsType::RealType RealType;
		
		MonteCarlo(const EngineParamsType& engineParams,
		           OperationsType& ops,
		           AlgorithmType& algorithm,
		           RandomNumberGeneratorType& rng) 
		 : engineParams_(engineParams),ops_(ops),rng_(rng),algorithm_(algorithm) { }

		PairType operator()(DynVarsType& dynVars,size_t iter)
		{
			PairType acc = PairType(0,0);
			ops_.set(dynVars);
			algorithm_.init();
			for (size_t j=0;j<dynVars.size;j++) {
				size_t i = ops_.proposeSite(j,rng_);
				ops_.proposeChange(i,rng_);
				bool flag= algorithm_.isAccepted(i,rng_,ops_,n);
				if (flag && !dynVars.isFrozen) { // Accepted
					algorithm_.accept(i,ops_);
					acc.first++;
				} else { // not accepted
				}
				acc.second++;
			} // lattice sweep
			return acc;
		}
		
	private:
		
		const EngineParamsType& engineParams_;
		OperationsType& ops_;
		RandomNumberGeneratorType& rng_;
		AlgorithmType& algorithm_;
		
	}; // MonteCarlo

	template<typename RngType,
	         typename ParametersType,
	         typename ModelType,
	         typename AlgorithmFactoryType,
	         typename ListType,
	         int n>
	class MonteCarloLoop
	{
		typedef std::pair<size_t,size_t> PairType;
		typedef typename Loki::TL::TypeAt<ListType,n>::Result OperationsType;
		typedef typename OperationsType::DynVarsType Type0;
		typedef MonteCarlo<ParametersType,ListType,AlgorithmFactoryType,
		                   RngType,n> MonteCarloType0;
		typedef typename ModelType::DynVarsType DynVarsType;

	public:

		static void loop(RngType& rng,
		                 const ParametersType& params,
		                 AlgorithmFactoryType& algorithm,
		                 ModelType& model,
		                 DynVarsType& dynVars,
		                 PsimagLite::Vector<PairType>::Type& accepted,
		                 size_t iter)
		{

			OperationsType* op = 0;
			model.setOperation(&op,n);
			MonteCarloType0 monteCarlo0(params,*op,algorithm,rng);
			Type0* spinPart = 0;
			dynVars.getField(&spinPart,n);
			PairType res= monteCarlo0(*spinPart,iter);
			accepted[n].first += res.first;
			accepted[n].second += res.second;
			MonteCarloLoop<RngType,ParametersType,ModelType,AlgorithmFactoryType,ListType,n-1>
			        ::loop(rng,params,algorithm,model,dynVars,accepted,iter);
		}
	};

	template<typename RngType,
	         typename ParametersType,
	         typename ModelType,
	         typename AlgorithmFactoryType,
	         typename ListType>
	class MonteCarloLoop<RngType,ParametersType,ModelType,AlgorithmFactoryType,ListType,-1>
	{
		typedef std::pair<size_t,size_t> PairType;
		typedef typename ModelType::DynVarsType DynVarsType;

	public:

		static void loop(RngType& rng,
		                 const ParametersType& params,
		                 AlgorithmFactoryType& algorithm,
		                 ModelType& model,
		                 DynVarsType& dynVars,
		                 PsimagLite::Vector<PairType>::Type& accepted,
		                 size_t iter)
		{}
	};

} // namespace Spf

/*@}*/
#endif // SPF_ENGINE_H

