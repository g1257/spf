/* DO NOT EDIT!!! Changes will be lost. Modify driver.pl instead
 This main.cpp was written by driver.pl
 SPF v7 by G.A. and C.S.
 Platform: linux
 Model: PnictidesTwoOrbitals
 */
#include "SimpleReader.h"
#include "ParametersEngine.h"
#include "Engine.h"
#include "ConcurrencySerial.h"
#include "ParametersPnictidesTwoOrbitals.h"
#include "PnictidesTwoOrbitals.h"
#include "GeometrySquare.h"
#include "RandomNumberGenerator.h"
#include "AlgorithmDiag.h"
#include "GreenFunction.h"

typedef double FieldType;
typedef Spf::ParametersEngine<FieldType> ParametersEngineType;
typedef Dmrg::ConcurrencySerial<FieldType> ConcurrencyType;
typedef Spf::GeometrySquare<FieldType> GeometryType;
typedef Spf::ParametersPnictidesTwoOrbitals<FieldType> ParametersModelType;
typedef Spf::PnictidesTwoOrbitals<ParametersEngineType,ParametersModelType,GeometryType> ModelType;
typedef ModelType::DynVarsType DynVarsType;
typedef Spf::RandomNumberGenerator<FieldType> RandomNumberGeneratorType;
typedef Spf::AlgorithmDiag<ParametersEngineType,ModelType,RandomNumberGeneratorType> AlgorithmType;
typedef Spf::GreenFunction<AlgorithmType> GreenFunctionType;
typedef Spf::Engine<ParametersEngineType,AlgorithmType,ModelType,ConcurrencyType,RandomNumberGeneratorType,GreenFunctionType> EngineType;

 
int main(int argc,char *argv[])
{
	ConcurrencyType concurrency(argc,argv);
	ParametersEngineType engineParams;
	ParametersModelType mp;
	Dmrg::SimpleReader reader(argv[1]);
	reader.load(engineParams);
	reader.load(mp);
	// print license
	std::string license = "Copyright © 2009 , UT-Battelle, LLC\n"
"All rights reserved\n"
"\n"
"[SPF, Version 7.0.0]\n"
"\n"
"*********************************************************\n"
"THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND\n"
"CONTRIBUTORS \"AS IS\" AND ANY EXPRESS OR IMPLIED\n"
"WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED\n"
"WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A\n"
"PARTICULAR PURPOSE ARE DISCLAIMED. \n"
"\n"
"Please see full open source license included in file LICENSE.\n"
"*********************************************************\n"
"\n"
"\n"
"// END LICENSE BLOCK\n"
;
	if (concurrency.root()) std::cerr<<license;
	GeometryType geometry(mp.linSize);
	
	ModelType model(engineParams,mp,geometry);
	AlgorithmType algorithm(engineParams,model);
	EngineType engine(engineParams,model,algorithm,concurrency);
	
	engine.main();
}


/* ####### End of main.cpp ######## */

