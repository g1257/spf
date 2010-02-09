// This file will be produced by the driver.pl later

#include "SimpleReader.h"
#include "ParametersEngine.h"
#include "Engine.h"
#include "ConcurrencySerial.h"
#include "ParametersPnictidesTwoOrbitals.h"
#include "PnictidesTwoOrbitals.h"
#include "GeometrySquare.h"

typedef double FieldType;
typedef Spf::ParametersEngine<FieldType> ParametersEngineType;
typedef Dmrg::ConcurrencySerial<FieldType> ConcurrencyType;
typedef Spf::GeometrySquare<FieldType> GeometryType;
typedef Spf::ParametersPnictidesTwoOrbitals<FieldType> ParametersModelType;
typedef Spf::PnictidesTwoOrbitals<FieldType,ParametersEngineType,ParametersModelType,GeometryType> ModelType;
typedef ModelType::DynVarsType DynVarsType;
typedef Spf::Engine<ParametersEngineType,ModelType,ConcurrencyType> EngineType;

 
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
	DynVarsType dynVars(geometry.volume(),engineParams.dynvarsfile);
	
	ModelType model(engineParams,mp,geometry);
	EngineType engine(engineParams,model,dynVars,concurrency);
	
	engine.main();
}
