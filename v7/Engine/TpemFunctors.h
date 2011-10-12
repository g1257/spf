
/** \ingroup SPF */
/*@{*/

/*! \file TpemFunctors.h
 *
 *  
 *
 */
#ifndef TPEM_FUNCTORS_H
#define TPEM_FUNCTORS_H

namespace Spf {
	template<typename TpemParametersType>
	class BaseFunctor {
	public:
		typedef TpemParametersType::RealType RealType;
		virtual RealType operator()(RealType x) const = 0;
	};
	
	template<typename TpemParametersType>
	class ActionFunctor : public BaseFunctor<RealType> {
	public:
		ActionFunctor(const TpemParametersType& tpemParameters)
		: tpemParameters_(tpemParameters)
		{}

		virtual RealType operator()(RealType x) const
		{
// 			tmpValues(a,b,mu,beta,1);
			return (a_ * x + b_ >= mu_) ?
			               log (1.0 + exp (-beta_ * (a_ * x + b_ - mu_)))
						   : -beta_ * (a_ * x + b_ - mu_)
				             + log (1.0 + exp (beta_ * (a_ * x + b_ - mu_)));
		}

	private:
		const TpemParametersType& tpemParameters_;
	}; // class ActionFunctor
	
	template<typename RealType>
	class EnergyFunctor  : public BaseFunctor<RealType> {
	public:
		EnergyFunctor(const TpemParametersType& tpemParameters)
		: tpemParameters_(tpemParameters)
		{}

		virtual RealType operator()(RealType x) const
		{
// 			tmpValues(a,b,mu,beta,1);
			return (a_ * x + b_) * 0.5 * 
			      (1.0 - tanh (0.5 * beta_ * (a_ * x + b_ - mu_)));
		}

	private:
		const TpemParametersType& tpemParameters_;
	}; // class EnergyFunctor
	
	template<typename RealType>
	class NumberFunctor : public BaseFunctor<RealType> {
	public:
		NumberFunctor(const TpemParametersType& tpemParameters)
		: tpemParameters_(tpemParameters)
		{}

		virtual RealType operator()(RealType x) const
		{
// 			tmpValues(a,b,mu,beta,1);
			return 0.5 * (1.0 - tanh (0.5 * beta_ * (a_ * x + b_ - mu_)));
		}
		
	private:
		const TpemParametersType& tpemParameters_;
	}; // class NumberFunctor
} // namespace Spf

/*@}*/
#endif // TPEM_FUNCTORS_H
