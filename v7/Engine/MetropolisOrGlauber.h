
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
		
		MetropolisOrGlauber(const PsimagLite::String& detailedBalance)
		{
			if (detailedBalance=="glauber") {
				algo_=GLAUBER;
			} else if (detailedBalance=="metropolis") {
				algo_=METROPOLIS;
			} else {
				PsimagLite::String s(__FILE__);
				s += " : Unknown detailed balance method: " + detailedBalance;
				s += "\n";
				throw PsimagLite::RuntimeError(s.c_str());
			}
		}
		
		bool operator()(const RealType& X2,RngType& rng) const
		{
			RealType X = X2;
			RealType r = rng();
			if (algo_==GLAUBER) {
				if (X<1) {
					X=X/(1.0+X);
				} else {
					X=1.0/(1.0+1.0/X);
				}
				return (X>r);
			}
			// METROPOLIS PROPER
			return (X > 1 || r < X);
		}

	private:
		bool algo_;
	}; // class MetropolisOrGlauber

} // namespace Spf

/*@}*/
#endif // METROPOLIS_OR_GLAUBER_H
