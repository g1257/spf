
/** \ingroup SPF */
/*@{*/

/*! \file MetropolisOrGlauber.h
 *
 *  
 *
 */
#ifndef METROPOLIS_OR_GLAUBER_H
#define METROPOLIS_OR_GLAUBER_H
#include <stdexcept>

namespace Spf {
	template<typename RealType,typename RngType>
	class MetropolisOrGlauber {
	public:
		enum { GLAUBER, METROPOLIS};
		
		MetropolisOrGlauber(const std::string& detailedBalance)
		{
			if (detailedBalance=="glauber") {
				algo_=GLAUBER;
			} else if (detailedBalance=="metropolis") {
				algo_=METROPOLIS;
			} else {
				std::string s(__FILE__);
				s += " : Unknown detailed balance method: " + detailedBalance;
				s += "\n";
				throw std::runtime_error(s.c_str());
			}
		}
		
		bool operator()(const RealType& X2,RngType& rng) const
		{
			RealType X = X2;
			if (algo_==GLAUBER) {
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
