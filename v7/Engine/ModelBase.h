
/** \ingroup SPF */
/*@{*/

/*! \file ModelBase.h
 *
 *  This is the interface that every model needs to implement
 *
 */
#ifndef MODEL_BASE_H
#define MODEL_BASE_H
#include <iostream>
#include <vector>
#include <cassert>

namespace Spf {
	template<typename DynVarsType,
	         typename EngineParamsType,
	         typename GeometryType>
	class ModelBase {
		typedef typename EngineParamsType::RealType RealType;

	public:

		// please follow the interface from one of the model classes
		// this isn't virtual due to performance reasons

	protected:

		void setTpemAandB(RealType& a,RealType &b,
		                  RealType& eMin,RealType& eMax) const
		{
			RealType factor = 1.02;
			if (eMin>0) eMin = 0;
			else eMin *= factor;
			if (eMax<0) throw std::runtime_error("Hmmm\n");
			else eMax *= factor;

			a = 0.5*(eMax-eMin);
			b = 0.5*(eMax+eMin);
// 			std::cerr<<"Set a="<<a<<" b="<<b<<"\n";
		}

		template<typename MatrixType>
		void setTpemSupport(PsimagLite::Vector<SizeType>::Type& support,
		                    const MatrixType& matrix,
		                    const MatrixType& matrix2,
		                    SizeType site) const
		{
			support.clear();
			for (SizeType i=0;i<matrix.n_col();i++) {
				typename MatrixType::value_type tmp = matrix(site,i) - matrix2(site,i);
				if (std::norm(tmp)<1e-6) continue;
				//if (std::find(support.begin(),support.end(),i)==
				//    support.end())
				support.push_back(i);
			}
			assert(support.size()>0);
// 			std::cerr<<"Set support with size="<<support.size()<<"\n";
// 			std::cerr<<support;
		}
	}; // ModelBase
} // namespace Spf

/*@}*/
#endif // MODEL_BASE_H
