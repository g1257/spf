%module myspf
%{
	#include "GeometrySquare.h"
	#include "Engine.h"
	#include "ConcurrencySerial.h"
	#include "ParametersEngine.h"
	#include "ParametersPnictidesThreeOrbitals.h"
	#include "PnictidesMultiOrbitals.h"
	#include "IoSimple.h"
	#include "Random48.h"
	#include "AlgorithmDiag.h"
	#include "GreenFunction.h"
%}	


namespace PsimagLite {
	template<typename T>
	class  Random48 {
	public:
		void seed(int long seed);
		T random();

					};//Class Random48
%template(Random48T) Random48<double>;
};// namespace Psimaglite


namespace Dmrg {
template<typename FieldType>
class ConcurrencySerial 
 
{
	public:
		
	ConcurrencySerial(int argc,char *argv[]);
		
	int nprocs();
		
	int rank();

	bool root();
	};//Class ConcurrencySerial
	
%template(ConcurrencyT) ConcurrencySerial<double>;

} // namespace Dmrg

namespace Spf {

template<typename T>
class GeometrySquare {
	public:
	GeometrySquare(SizeType l);

	SizeType length() const;

					};//Class GeometrySquare

%template(GeometryT) GeometrySquare<double>;

 
class IoSimpleIn {
	public:
		IoSimpleIn(const char *fn);
		const char* filename() const;
	
				};//Class IoSimpleIn


template<typename FieldType_,typename IoInType>
	struct ParametersEngine {
		
		//! Read Dmrg parameters from inp file
		ParametersEngine(IoInType& io);
		
		SizeType latticeLength;
	
							};//Struct ParametersEngine

%template(ParametersEngineT) ParametersEngine<double,IoSimpleIn>;


template<typename ParametersEngineType,typename IoInType>
	struct ParametersPnictidesThreeOrbitals {
		

		ParametersPnictidesThreeOrbitals(IoInType& io,
                                const ParametersEngineType& engineParams);
	
											};//Struct ParametersPnictidesThreeOrbitals
%template(ParametersModelT) ParametersPnictidesThreeOrbitals<ParametersEngine<double,IoSimpleIn>,IoSimpleIn>;


template<typename EngineParamsType,typename ParametersModelType,typename GeometryType,typename ConcurrencyType>
	class PnictidesMultiOrbitals { 
		public:
 
	PnictidesMultiOrbitals(const EngineParamsType& engineParams,
				const ParametersModelType& mp,
				const GeometryType& geometry,
                                ConcurrencyType& concurrency);

								};//Class PnictidesMultiOrbitals
%template(ModelT) PnictidesMultiOrbitals<
	ParametersEngine<double,IoSimpleIn>,
	ParametersPnictidesThreeOrbitals<ParametersEngine<double,IoSimpleIn>,IoSimpleIn>,
	GeometrySquare<double>,
	Dmrg::ConcurrencySerial<double>
>;


template<typename EngineParametersType,typename ModelType,typename RngType>
	class AlgorithmDiag {
	public:
			AlgorithmDiag(const EngineParametersType& engineParams,ModelType& model);
	

						};//Class AlgorithmDiag 
%template(AlgorithmT) AlgorithmDiag<ParametersEngine<double,IoSimpleIn>,
PnictidesMultiOrbitals<
	ParametersEngine<double,IoSimpleIn>,
	ParametersPnictidesThreeOrbitals<ParametersEngine<double,IoSimpleIn>,IoSimpleIn>,
	GeometrySquare<double>,Dmrg::ConcurrencySerial<double> >,
PsimagLite::Random48<double> >;


template<typename EngineParametersType,typename AlgorithmType>
	class GreenFunction {

	public:
		typedef typename AlgorithmType::FieldType FieldType;
		typedef typename AlgorithmType::ComplexType ComplexType;

		GreenFunction(const EngineParametersType& engineParams,AlgorithmType& algorithm,SizeType hilbertSize) :
			engineParams_(engineParams),algorithm_(algorithm),hilbertSize_(hilbertSize),data_(hilbertSize,hilbertSize);

						};//Class GreenFunction

%template(GreenFunctionD) GreenFunction<ParametersEngine<double,IoSimpleIn>,
AlgorithmDiag<ParametersEngine<double,IoSimpleIn>,
	PnictidesMultiOrbitals<ParametersEngine<double,IoSimpleIn>,
		ParametersPnictidesThreeOrbitals<ParametersEngine<double,IoSimpleIn>,IoSimpleIn>,
	GeometrySquare<double>,Dmrg::ConcurrencySerial<double> >,
	PsimagLite::Random48<double> > >;



template<typename ParametersType,typename AlgorithmType,typename ModelType,
 		typename ConcurrencyType,typename RandomNumberGeneratorType,typename GreenFunctionType>
	class Engine {
		
		typedef typename ParametersType::FieldType FieldType;
		typedef typename ModelType::DynVarsType DynVarsType;
		typedef Dmrg::ProgressIndicator ProgressIndicatorType;
		typedef std::pair<SizeType,SizeType> PairType;
		
		public:
			
		Engine(ParametersType& params,ModelType& model,AlgorithmType& algorithm,ConcurrencyType& concurrency);

void main();
				}; //Class Engine
%template(EngineT) Engine<
ParametersEngine<double,IoSimpleIn>,
AlgorithmDiag<ParametersEngine<double,IoSimpleIn>,
	PnictidesMultiOrbitals<
	ParametersEngine<double,IoSimpleIn>,
	ParametersPnictidesThreeOrbitals<ParametersEngine<double,IoSimpleIn>,IoSimpleIn>,
	GeometrySquare<double>,Dmrg::ConcurrencySerial<double> >,
	PsimagLite::Random48<double> >,
PnictidesMultiOrbitals<
	ParametersEngine<double,IoSimpleIn>,
	ParametersPnictidesThreeOrbitals<ParametersEngine<double,IoSimpleIn>,IoSimpleIn>,
	GeometrySquare<double>,Dmrg::ConcurrencySerial<double> >,
Dmrg::ConcurrencySerial<double>,
PsimagLite::Random48<double>,
GreenFunction<ParametersEngine<double,IoSimpleIn>,
AlgorithmDiag<ParametersEngine<double,IoSimpleIn>,
	PnictidesMultiOrbitals<ParametersEngine<double,IoSimpleIn>,
	ParametersPnictidesThreeOrbitals<ParametersEngine<double,IoSimpleIn>,IoSimpleIn>,
	GeometrySquare<double>,Dmrg::ConcurrencySerial<double> >,
	PsimagLite::Random48<double> > >
>;

}//namespace Spf

