// BEGIN LICENSE BLOCK
/*
Copyright (c) 2009 , UT-Battelle, LLC
All rights reserved

[DMRG++, Version 2.0.0]
[by G.A., Oak Ridge National Laboratory]

UT Battelle Open Source Software License 11242008

OPEN SOURCE LICENSE

Subject to the conditions of this License, each
contributor to this software hereby grants, free of
charge, to any person obtaining a copy of this software
and associated documentation files (the "Software"), a
perpetual, worldwide, non-exclusive, no-charge,
royalty-free, irrevocable copyright license to use, copy,
modify, merge, publish, distribute, and/or sublicense
copies of the Software.

1. Redistributions of Software must retain the above
copyright and license notices, this list of conditions,
and the following disclaimer.  Changes or modifications
to, or derivative works of, the Software should be noted
with comments and the contributor and organization's
name.

2. Neither the names of UT-Battelle, LLC or the
Department of Energy nor the names of the Software
contributors may be used to endorse or promote products
derived from this software without specific prior written
permission of UT-Battelle.

3. The software and the end-user documentation included
with the redistribution, with or without modification,
must include the following acknowledgment:

"This product includes software produced by UT-Battelle,
LLC under Contract No. DE-AC05-00OR22725  with the
Department of Energy."
 
*********************************************************
DISCLAIMER

THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT OWNER, CONTRIBUTORS, UNITED STATES GOVERNMENT,
OR THE UNITED STATES DEPARTMENT OF ENERGY BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.

NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED
STATES DEPARTMENT OF ENERGY, NOR THE COPYRIGHT OWNER, NOR
ANY OF THEIR EMPLOYEES, REPRESENTS THAT THE USE OF ANY
INFORMATION, DATA, APPARATUS, PRODUCT, OR PROCESS
DISCLOSED WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*********************************************************


*/
// END LICENSE BLOCK
/** \ingroup SPF */
/*@{*/

/*! \file ParametersEngine.h
 *
 *  Contains the parameters for the DmrgSolver class and implements functionality to
 *  read them from a JSON file
 *
 */
#ifndef PARAMETERSENGINE_H
#define PARAMETERSENGINE_H
#include <string>
#include <vector>
#include <iostream>
#include <cstdlib>
#include "TypeToString.h"

namespace Spf {
	
	//! Structure that contains the Engine parameters
	template<typename RealType_,typename IoInType_>
	struct ParametersEngine {
		typedef RealType_ RealType;
		typedef IoInType_ IoInType;
		//! Read Dmrg parameters from inp file
		ParametersEngine(IoInType& io)
		{
			io.readline(options,"EngineOptions=");
			io.readline(version,"Version=");
			io.readline(filename,"OutputFilename=");
			io.readline(carriers,"Carriers=");
			io.readline(mu,"ChemicalPotential=");
			io.readline(beta,"Beta=");
			io.readline(iterTherm,"MonteCarloThermalizations=");
			io.readline(iterEffective,"MonteCarloEffectiveIteractions=");
			io.readline(iterUnmeasured,"MonteCarloUnmeasuredPlusOne=");
			io.read(mcWindow,"MonteCarloWindows");
			io.readline(dynvarsfile,"MonteCarloStartTypeOrFile=");
			io.readline(dynvarslevel,"MonteCarloStartLevel=");
			io.readline(histSteps,"HistogramSteps=");
			io.readline(boundaryConditions,"BoundaryConditions=");
			std::string s;
			io.readline(s,"RandomSeed=");
			if (s == "TIME" || s == "time")
				randomSeed = -1;
			else    randomSeed = atoi(s.c_str());
			io.readline(latticeLength,"LatticeLength=");
			
			coresForKernel = 1;
			try {
				io.readline(coresForKernel,"CoresForKernel=");
			} catch (std::exception& e) {
				io.rewind();
			}
			saveEach=0;
			try {
				io.readline(saveEach,"SaveEach=");
			} catch (std::exception& e) {
				io.rewind();
			}
		}

		std::string filename; // filename to save observables and continued fractions
		std::string version;
		std::string options; // options
		size_t carriers;
		mutable RealType mu; // chemical potential 
		RealType beta; // inverse temperature
		size_t iterTherm,iterEffective,iterUnmeasured;
		std::vector<RealType> mcWindow; // monte carlo window of change
		std::string dynvarsfile; // file with fields to start from or none
		size_t dynvarslevel; // from where to start to read in dynvarsfile
		size_t histSteps; // histogram steps
		std::string boundaryConditions; // boundary conditions
		long int randomSeed;
		size_t latticeLength;
		size_t coresForKernel;
		size_t saveEach;
	};

	//! print dmrg parameters
	template<typename RealType,typename IoInType>
	std::ostream &operator<<(std::ostream &os,
		ParametersEngine<RealType,IoInType> const &parameters)
	{
		os<<"#This is SPF\n";
		os<<"parameters.version="<<parameters.version<<"\n";
		os<<"parameters.filename="<<parameters.filename<<"\n";
		os<<"parameters.options="<<parameters.options<<"\n";
		os<<"parameters.carriers="<<parameters.carriers<<"\n";
		os<<"parameters.mu="<<parameters.mu<<"\n";
		os<<"parameters.beta="<<parameters.beta<<"\n";
		os<<"parameters.iterTherm="<<parameters.iterTherm<<"\n";
		os<<"parameters.iterEffective="<<parameters.iterEffective<<"\n";
		os<<"parameters.iterUnmeasured="<<parameters.iterUnmeasured<<"\n";
		os<<"parameters.mcWindow\n";
		os<<parameters.mcWindow;
		os<<"parameters.dynvarsfile="<<parameters.dynvarsfile<<"\n";
		os<<"parameters.dynvarslevel="<<parameters.dynvarslevel<<"\n";
		os<<"parameters.histSteps="<<parameters.histSteps<<"\n";
		os<<"parameters.boundaryConditions="<<parameters.boundaryConditions<<"\n";
		std::string s = ttos(parameters.randomSeed);
		if (parameters.randomSeed<0) s="TIME";
		os<<"parameters.randomSeed="<<s<<"\n";
		os<<"parameters.latticeLength="<<parameters.latticeLength<<"\n";
		os<<"parameters.coresForKernel="<<parameters.coresForKernel<<"\n";
		return os;
	}
} // namespace Dmrg
/*@}*/

#endif
