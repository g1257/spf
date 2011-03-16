
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
	template<typename DynVarsType,typename EngineParamsType,typename ParametersModelType_,typename GeometryType>
	class ModelBase {
		typedef typename EngineParamsType::FieldType FieldType;
		
		public:
		
		size_t totalFlips() const;
		
		size_t hilbertSize() const;
		
		FieldType deltaDirect(size_t i) const;
	
		void set(DynVarsType& dynVars);
		
		template<typename RandomNumberGeneratorType>
		void propose(size_t i,RandomNumberGeneratorType& rng);
				
		void doMeasurements(DynVarsType& dynVars, size_t iter,std::ostream& fout);
		
		void fillAndDiag(std::vector<FieldType> &eig);
		
		void fillAndDiag(std::vector<FieldType> &eig,const DynVarsType& dynVars,char jobz='N');
		
		void adjustChemPot(const std::vector<FieldType>& eigs);
		
		void accept(size_t i);
		
		FieldType integrationMeasure(size_t i);
	}; // ModelBase
} // namespace Spf

/*@}*/
#endif // MODEL_BASE_H
