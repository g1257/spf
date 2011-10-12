
/** \ingroup TPEM */
/*@{*/

/*! \file TpemFunctors.h
 *
 *  
 *
 */
#ifndef TPEM_FUNCTORS_H
#define TPEM_FUNCTORS_H
#include <cmath>

namespace Tpem {
	template<typename TpemParametersType>
	class BaseFunctor {
	public:
		typedef typename TpemParametersType::RealType RealType;
		virtual RealType operator()(RealType x) const = 0;
	};
	
	template<typename TpemParametersType>
	class ActionFunctor : public BaseFunctor<TpemParametersType> {
	public:
		typedef typename TpemParametersType::RealType RealType;

		ActionFunctor(const TpemParametersType& tpemParameters)
		: tpemParameters_(tpemParameters)
		{}

		virtual RealType operator()(RealType x) const
		{
// 			tmpValues(a,b,mu,beta,1);
			RealType tmp = tpemParameters_.a * x + tpemParameters_.b - tpemParameters_.mu;
			const RealType& beta = tpemParameters_.beta;
			return (tmp>=0) ? log (1.0 + exp (beta*tmp)) 
			                : -beta * tmp + log (1.0 + exp (beta * tmp));
		}

	private:
		const TpemParametersType& tpemParameters_;
	}; // class ActionFunctor
	
	template<typename TpemParametersType>
	class EnergyFunctor  : public BaseFunctor<TpemParametersType> {
	public:
		typedef typename TpemParametersType::RealType RealType;

		EnergyFunctor(const TpemParametersType& tpemParameters)
		: tpemParameters_(tpemParameters)
		{}

		virtual RealType operator()(RealType x) const
		{
// 			tmpValues(a,b,mu,beta,1);
			RealType tmp = tpemParameters_.a * x + tpemParameters_.b;
			return tmp*0.5*(1.0-tanh(0.5*tpemParameters_.beta*(tmp-tpemParameters_.mu)));
		}

	private:
		const TpemParametersType& tpemParameters_;
	}; // class EnergyFunctor
	
	template<typename TpemParametersType>
	class NumberFunctor : public BaseFunctor<TpemParametersType> {
	public:
		typedef typename TpemParametersType::RealType RealType;

		NumberFunctor(const TpemParametersType& tpemParameters)
		: tpemParameters_(tpemParameters)
		{}

		virtual RealType operator()(RealType x) const
		{
// 			tmpValues(a,b,mu,beta,1);
			RealType tmp = tpemParameters_.a * x + tpemParameters_.b - tpemParameters_.mu;
			return 0.5*(1.0-tanh(0.5*tpemParameters_.beta*tmp));
		}
		
	private:
		const TpemParametersType& tpemParameters_;
	}; // class NumberFunctor
} // namespace Tpem

/*@}*/
#endif // TPEM_FUNCTORS_H
