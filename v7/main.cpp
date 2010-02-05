// This file will be produced by the driver.pl later

#include "SimpleReader.h"
#include "ParametersEngine.h"
#include "ConcurrencySerial.h"
#include "ParametersPnictidesTwoOrbitals.h"

typedef double FieldType;
typedef Dmrg::ConcurrencySerial<FieldType> ConcurrencyType;
typedef Spf::ParametersPnictidesTwoOrbitals<FieldType> ParametersModelType;
 
int main(int argc,char *argv[])
{
	ConcurrencyType concurrency(argc,argv);
	ParametersModelType mp;
	Spf::ParametersEngine<FieldType> engineParams;
	Dmrg::SimpleReader reader(argv[1]);
	reader.load(mp);
	reader.load(engineParams);
	// print license
	std::string license = "Copyright © 2009 , UT-Battelle, LLC\n"
"All rights reserved\n"
"\n"
"[DMRG++, Version 2.0.0]\n"
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
	/*GeometryType geometry;
	ModelType model(mp,geometry);
	EngineType engine(engineParams,model,concurrency);
	
	engine.main(); */
}
