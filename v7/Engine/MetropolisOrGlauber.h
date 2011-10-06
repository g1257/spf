
/** \ingroup SPF */
/*@{*/

/*! \file MetropolisOrGlauber.h
 *
 *  
 *
 */
#ifndef METROPOLIS_OR_GLAUBER_H
#define METROPOLIS_OR_GLAUBER_H

namespace Spf {
	template<typename RealType,typename RngType>
	class MetropolisOrGlauber {
	public:
		enum { GLAUBER, METROPOLIS};
		
		MetropolisOrGlauber(size_t which = GLAUBER) : algo_(which)
		{}
		
		bool operator()(const RealType& X2,RngType& rng) const
		{
			RealType X = X2;
			if (algo_) {
				if (X<1) {
					X=X/(1.0+X);
				} else {
					X=1.0/(1.0+1.0/X);
				}
				return (X>rng());
			}
			// METROPOLIS PROPER
			return (X > 1 || rng() < X);
		}

	private:
		bool algo_;
	}; // class MetropolisOrGlauber

} // namespace Spf

/*@}*/
#endif // METROPOLIS_OR_GLAUBER_H
