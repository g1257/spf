
/** \ingroup TPEM */
/*@{*/

/*! \file TpemParameters.h
 *
 *  Do not add functions here except a ctor
 *
 */
#ifndef TPEM_PARAMETERS_H
#define TPEM_PARAMETERS_H
#include "CrsMatrix.h"
#include <cassert>

namespace Tpem {
	
	template<typename RealType>
	class EmptyCallBack {
	public:
		void setTpemThings(RealType& a, RealType& b, std::vector<size_t>& support) const
		{}
	};
	
	template<typename IoInType,typename RealType_,typename CallbackType=EmptyCallBack<RealType_> >
	struct TpemParameters {
		typedef RealType_ RealType;
		enum {TPEM,PEM};
		
		TpemParameters(IoInType& io,
		               const RealType& mu1,
		               const RealType& beta1,
		               const CallbackType* callback=0)
		: mu(mu1),beta(beta1)
		{
			std::string s;
			io.readline(s,"TpemAlgorithm=");
			if (s=="TPEM" or s=="tpem" or s=="Tpem") {
				algorithm=TPEM;
			} else if (s=="PEM" or s=="pem" or s=="Pem") {
				algorithm=PEM;
			} else {
				s = std::string("Expected tpem or pem but found ") + s;
				throw std::runtime_error(s.c_str());
			}
			io.readline(cutoff,"TpemCutoff=");
			io.readline(epsForProduct,"TpemEpsForProduct=");
			io.readline(epsForTrace,"TpemEpsForTrace=");
			if (callback) {
				callback->setTpemThings(a,b,support);
				return;
			}
			io.readline(a,"TpemA=");
			io.readline(b,"TpemB=");
			io.read(support,"TpemSupport");
		}

		const RealType& mu;
		const RealType& beta;
		size_t cutoff;
		size_t algorithm;
		RealType epsForProduct;
		RealType epsForTrace;
		RealType a;
		RealType b;
		std::vector<size_t> support;
	}; // struct TpemParameters
	
	template<typename IoInType,typename RealType,typename CallbackType>
	std::ostream& operator<<(std::ostream& os,const TpemParameters<IoInType,RealType,CallbackType>& t)
	{
		os<<"#TpemParameters\n";
		os<<"#TpemMu="<<t.mu<<"\n";
		os<<"#TpemBeta="<<t.beta<<"\n";
		os<<"#TpemCutoff="<<t.cutoff<<"\n";
		os<<"#TpemAlgorithm="<<t.algorithm<<"\n";
		os<<"#TpemEpsProduct="<<t.epsForProduct<<"\n";
		os<<"#TpemEpsTrace="<<t.epsForTrace<<"\n";
		os<<"#TpemA="<<t.a<<"\n";
		os<<"#TpemB="<<t.b<<"\n";
		os<<"#TpemSupport\n";
		os<<t.support;
		return os;
	}
} // namespace Tpem

/*@}*/
#endif // TPEM_PARAMETERS_H
