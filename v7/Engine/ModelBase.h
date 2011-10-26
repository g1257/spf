
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

namespace Spf {
	template<typename DynVarsType,
	         typename EngineParamsType,
	         typename ParametersModelType_,
	         typename GeometryType,
	         typename ConcurrencyType>
	class ModelBase {
		typedef typename EngineParamsType::RealType RealType;
		
	public:

		size_t totalFlips() const;

		size_t hilbertSize() const;

		RealType deltaDirect(size_t i) const;
	
		void set(DynVarsType& dynVars);

		template<typename RandomNumberGeneratorType>
		void propose(size_t i,RandomNumberGeneratorType& rng);

		void doMeasurements(DynVarsType& dynVars, size_t iter,std::ostream& fout);

		void fillAndDiag(std::vector<RealType> &eig);

		void fillAndDiag(std::vector<RealType> &eig,const DynVarsType& dynVars,char jobz='N');

		void adjustChemPot(const std::vector<RealType>& eigs);

		void accept(size_t i);

		RealType integrationMeasure(size_t i);

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
			std::cerr<<"Set a="<<a<<" b="<<b<<"\n";
		}

		template<typename MatrixType>
		void setTpemSupport(std::vector<size_t>& support,
		                    const MatrixType& matrix,
		                    const MatrixType& matrix2,
		                    size_t site) const
		{
			support.clear();
			for (size_t i=0;i<matrix.n_col();i++) {
				typename MatrixType::value_type tmp = matrix(site,i) - matrix2(site,i);
				if (std::norm(tmp)<1e-6) continue;
				//if (std::find(support.begin(),support.end(),i)==
				//    support.end())
				support.push_back(i);
			}
			assert(support.size()>0);
			std::cerr<<"Set support with size="<<support.size()<<"\n";
			std::cerr<<support;
		}

	}; // ModelBase
} // namespace Spf

/*@}*/
#endif // MODEL_BASE_H
