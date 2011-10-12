
/** \ingroup SPF */
/*@{*/

/*! \file ModelBase.h
 *
 *  This is the interface that every model needs to implement
 *
 */
#ifndef MODEL_BASE_H
#define MODEL_BASE_H


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
	}; // ModelBase
} // namespace Spf

/*@}*/
#endif // MODEL_BASE_H
