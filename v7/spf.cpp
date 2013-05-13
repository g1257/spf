#include "String.h"
#include <iostream>
void printLicense()
{
	const PsimagLite::String license_=
	        "Copyright (c) 2009-2013, UT-Battelle, LLC\n"
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
	        "\n";
	std::cout<<license_;
}

typedef double FieldType;
#include "ParametersEngine.h"
#include "Engine.h"
#include <algorithm>
#ifndef USE_MPI
#include "ConcurrencySerial.h"
typedef PsimagLite::ConcurrencySerial<FieldType> MyConcurrencyType;
#else
#include "ConcurrencyMpi.h"
typedef PsimagLite::ConcurrencyMpi<FieldType> MyConcurrencyType;
#endif
#include "PnictidesMultiOrbitals.h"
#include "DmsMultiOrbital.h"
#include "PhononsTwoOrbitals.h"
#include "GeometrySquare.h"
#include "GeometryCubic.h"
#include "GeometryFcc.h"
#include "GeometrySquare45Degrees.h"
#include "Random48.h"

typedef PsimagLite::IoSimple::In IoInType;
typedef Spf::ParametersEngine<FieldType,IoInType,MyConcurrencyType> ParametersEngineType;
typedef Spf::GeometrySquare<FieldType> GeometrySquareType;
typedef Spf::GeometryCubic<FieldType> GeometryCubicType;
typedef Spf::GeometryFcc<FieldType> GeometryFccType;
typedef Spf::GeometrySquare45Degrees<FieldType> GeometrySquare45DegreesType;

//typedef ModelType::DynVarsType DynVarsType;
typedef PsimagLite::Random48<FieldType> RandomNumberGeneratorType;

template<typename GeometryType,typename ModelType>
void mainLoop2(const ParametersEngineType& engineParams,
               IoInType& io,
               const GeometryType& geometry,
               MyConcurrencyType& concurrency)
{
	typedef Spf::Engine<ParametersEngineType,ModelType,IoInType,RandomNumberGeneratorType> EngineType;

	ModelType model(engineParams,io,geometry,concurrency);

	EngineType engine(engineParams,model,io,concurrency);

	engine.main();
}

template<typename GeometryType>
void mainLoop(const ParametersEngineType& engineParams,
              IoInType& io,
              MyConcurrencyType& concurrency)
{
	typedef Spf::PnictidesMultiOrbitals<ParametersEngineType,GeometryType,MyConcurrencyType> PnictidesMultiOrbitalsType;
	typedef Spf::DmsMultiOrbital<ParametersEngineType,GeometryType,MyConcurrencyType> DmsMultiOrbitalType;
	typedef Spf::PhononsTwoOrbitals<ParametersEngineType,GeometryType,MyConcurrencyType> PhononsTwoOrbitalsType;

	GeometryType geometry(engineParams.latticeLength);

	if (engineParams.model=="DmsMultiOrbital") {
		mainLoop2<GeometryType,DmsMultiOrbitalType>(engineParams,io,geometry,concurrency);
	} else if (engineParams.model=="PnictidesMultiOrbitals") {
		mainLoop2<GeometryType,PnictidesMultiOrbitalsType>(engineParams,io,geometry,concurrency);
	} else if (engineParams.model=="PhononsTwoOrbitals") {
		mainLoop2<GeometryType,PhononsTwoOrbitalsType>(engineParams,io,geometry,concurrency);
	} else {
		std::cerr<<"model="<<engineParams.model<<"\n";
		throw PsimagLite::RuntimeError("Unknown model");
	}
}

int main(int argc,char *argv[])
{
	MyConcurrencyType concurrency(argc,argv);
	if (argc<2) {
		PsimagLite::String s = "Usage is: ./" + PsimagLite::String(argv[0]) +
		" input_file\n";
		throw PsimagLite::RuntimeError(s.c_str());
	}
	PsimagLite::IoSimple::In io(argv[1]);
	ParametersEngineType engineParams(io,concurrency);
	
	if (concurrency.root()) printLicense();
	if (engineParams.geometry=="square") {
		mainLoop<GeometrySquareType>(engineParams,io,concurrency);
	} else if (engineParams.geometry=="cubic") {
		mainLoop<GeometryCubicType>(engineParams,io,concurrency);
	} else if (engineParams.geometry=="fcc") {
		mainLoop<GeometryFccType>(engineParams,io,concurrency);
	} else if (engineParams.geometry=="square45Degrees") {
		mainLoop<GeometrySquare45DegreesType>(engineParams,io,concurrency);
	} else {
		std::cerr<<"geometry="<<engineParams.geometry<<"\n";
		throw PsimagLite::RuntimeError("Unknown geometry\n");
	}
}


/* ####### End of spf.cpp ######## */

