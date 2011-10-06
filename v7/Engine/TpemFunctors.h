
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

	template<typename RealType>
	class ActionFunctor {
	public:
		ActionFunctor(const RealType& a,
					  const RealType& b,
				const RealType& mu,
				const RealType& beta)
		: a_(a),b_(b),mu_(mu),beta_(beta)
		{}

		RealType operator()(RealType x)
		{
// 			tmpValues(a,b,mu,beta,1);
			return (a_ * x + b_ >= mu_) ?
			               log (1.0 + exp (-beta_ * (a_ * x + b_ - mu_)))
						   : -beta_ * (a_ * x + b_ - mu_)
				             + log (1.0 + exp (beta_ * (a_ * x + b_ - mu_)));
		}

	private:
		const RealType& a_;
		const RealType& b_;
		const RealType& mu_;
		const RealType& beta_;
	}; // class ActionFunctor
	
	template<typename RealType>
	class EnergyFunctor  {
	public:
		EnergyFunctor(const RealType& a,
		               const RealType& b,
		               const RealType& mu,
		               const RealType& beta)
		: a_(a),b_(b),mu_(mu),beta_(beta)
		{}

		RealType operator()(RealType x)
		{
// 			tmpValues(a,b,mu,beta,1);
			return (a_ * x + b_) * 0.5 * 
			      (1.0 - tanh (0.5 * beta_ * (a_ * x + b_ - mu_)));
		}

	private:
		const RealType& a_;
		const RealType& b_;
		const RealType& mu_;
		const RealType& beta_;
	}; // class EnergyFunctor
	
	template<typename RealType>
	class NumberFunctor {
	public:
		NumberFunctor(const RealType& a,
					  const RealType& b,
				const RealType& mu,
				const RealType& beta)
		: a_(a),b_(b),mu_(mu),beta_(beta)
		{}

		RealType operator()(RealType x)
		{
// 			tmpValues(a,b,mu,beta,1);
			return 0.5 * (1.0 - tanh (0.5 * beta_ * (a_ * x + b_ - mu_)));
		}
		
	private:
		const RealType& a_;
		const RealType& b_;
		const RealType& mu_;
		const RealType& beta_;
	}; // class NumberFunctor
} // namespace Spf

/*@}*/
#endif // TPEM_FUNCTORS_H
