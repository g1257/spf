/* Copyright (c) 2009-2013, UT-Battelle, LLC
All rights reserved

[SPF, Version 7.0]
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
/** \ingroup SPF */
/*@{*/

/*! \file ParametersEngine.h
 *
 *
 */
#ifndef PARAMETERSENGINE_H
#define PARAMETERSENGINE_H
#include "String.h"
#include <vector>
#include <iostream>
#include <cstdlib>
#include <stdexcept>
#include "Map.h"
#include "TypeToString.h"
#include "IoSimple.h"
#include "Concurrency.h"

namespace Spf {

	//! Structure that contains the Engine parameters
	template<typename RealType_,typename IoInType_>
	struct ParametersEngine {
		typedef RealType_ RealType;
		typedef IoInType_ IoInType;
		ParametersEngine(IoInType& io)
		{
			io.readline(options,"EngineOptions=");
			io.readline(geometry,"GeometryKind=");
			io.readline(model,"Model=");
			io.readline(version,"Version=");
			io.readline(filename,"OutputFilename=");
			int temp = 0;
			io.readline(temp,"TargetElectronsUp=");
			carriers = temp;
			io.readline(temp,"TargetElectronsDown=");
			carriers += temp;
			io.readline(mu,"ChemicalPotential=");
			io.readline(beta,"Beta=");
			io.readline(iterTherm,"MonteCarloThermalizations=");
			io.readline(iterEffective,"MonteCarloEffectiveIteractions=");
			io.readline(iterUnmeasured,"MonteCarloUnmeasuredPlusOne=");

			bool flag = false;
			try {
				typename PsimagLite::Vector<RealType>::Type tmpVector;
				io.read(tmpVector,"MonteCarloWindows");
				flag = true;
			} catch (std::exception& e) {}

			if (flag) {
				PsimagLite::String str("MonteCarloWindows is not supported anymore");
				str += " Please add MonteCarloWindow[\"SpinPhi\"]= and ";
				str += " MonteCarloWindow[\"SpinTheta\"]= to the input file\n";
				throw PsimagLite::RuntimeError(str);
			}

			io.read(mcWindow,"MonteCarloWindow");
			io.readline(dynvarsfile,"MonteCarloStartTypeOrFile=");
			io.readline(dynvarslevel,"MonteCarloStartLevel=");
			io.readline(histSteps,"HistogramSteps=");
			io.readline(boundaryConditions,"BoundaryConditions=");
			PsimagLite::String s;
			io.readline(s,"RandomSeed=");
			if (s == "TIME" || s == "time")
				randomSeed = -1;
			else    randomSeed = atoi(s.c_str());
			io.readline(latticeLength,"LadderLeg=");

			coresForKernel = 1;
			try {
				io.readline(coresForKernel,"CoresForKernel=");
			} catch (std::exception& e) {

				if (PsimagLite::Concurrency::nprocs()>1) {
					std::cerr<<"Did you forget CoresForKernel= line in the input file?\n";
					throw e;
				}
			}
			if (SizeType(PsimagLite::Concurrency::nprocs())<coresForKernel) {
				s= PsimagLite::String(__FILE__) + " " + ttos(__LINE__) + " " + __FUNCTION__;
				s += PsimagLite::String(": nprocs<coresForKernel is an error\n");
				throw PsimagLite::RuntimeError(s.c_str());
			}
			saveEach=0;

//			io.rewind();
			try {
				io.readline(saveEach,"SaveEach=");
			} catch (std::exception& e) {
//				io.rewind();
			}
			detailedBalance = "glauber";
			try {
				io.readline(detailedBalance,"DetailedBalance=");
			} catch (std::exception& e) {
				PsimagLite::String s("*** WARNING *** Omission of DetailedBalance line ");
				s += " in the input file is deprecated. Assuming GLAUBER. ";
				s += "Please make sure it is correct.\n";
				std::cerr<<s;
			}
//			io.rewind();

			adjustEach = 0;
			try {
				io.readline(adjustEach,"AdjustEach=");
			} catch (std::exception& e) {
				PsimagLite::String s("*** WARNING *** Omission of AdjustEach= line ");
				s += " in the input file is deprecated. Assuming 0. ";
				s += "Please make sure it is correct.\n";
				std::cerr<<s;
			}
			npthreads = 1;
//			io.rewind();

		}

		PsimagLite::String filename; // filename to save observables and continued fractions
		PsimagLite::String geometry;
		PsimagLite::String model;
		PsimagLite::String version;
		PsimagLite::String options; // options
		RealType carriers;
		mutable RealType mu; // chemical potential
		RealType beta; // inverse temperature
		SizeType iterTherm,iterEffective,iterUnmeasured;
		mutable typename PsimagLite::Map<PsimagLite::String,RealType>::Type  mcWindow; // monte carlo windows
		PsimagLite::String dynvarsfile; // file with fields to start from or none
		SizeType dynvarslevel; // from where to start to read in dynvarsfile
		SizeType histSteps; // histogram steps
		PsimagLite::String boundaryConditions; // boundary conditions
		long int randomSeed;
		SizeType latticeLength;
		SizeType coresForKernel;
		SizeType saveEach;
		PsimagLite::String detailedBalance;
		SizeType adjustEach;
		SizeType npthreads;
	};

	//! print dmrg parameters
	template<typename RealType,typename IoInType>
	std::ostream &operator<<(std::ostream &os,
		ParametersEngine<RealType,IoInType> const &parameters)
	{
		os<<"#This is SPF\n";
		os<<"parameters.version="<<parameters.version<<"\n";
		os<<"parameters.geometry="<<parameters.geometry<<"\n";
		os<<"parameters.model="<<parameters.model<<"\n";
		os<<"parameters.filename="<<parameters.filename<<"\n";
		os<<"parameters.options="<<parameters.options<<"\n";
		os<<"parameters.carriers="<<parameters.carriers<<"\n";
		os<<"parameters.mu="<<parameters.mu<<"\n";
		os<<"parameters.beta="<<parameters.beta<<"\n";
		os<<"parameters.iterTherm="<<parameters.iterTherm<<"\n";
		os<<"parameters.iterEffective="<<parameters.iterEffective<<"\n";
		os<<"parameters.iterUnmeasured="<<parameters.iterUnmeasured<<"\n";
		PsimagLite::printMap(os,parameters.mcWindow,"parameters.mcWindow");
		os<<"parameters.dynvarsfile="<<parameters.dynvarsfile<<"\n";
		os<<"parameters.dynvarslevel="<<parameters.dynvarslevel<<"\n";
		os<<"parameters.histSteps="<<parameters.histSteps<<"\n";
		os<<"parameters.boundaryConditions="<<parameters.boundaryConditions<<"\n";
		PsimagLite::String s = ttos(parameters.randomSeed);
		if (parameters.randomSeed<0) s="TIME";
		os<<"parameters.randomSeed="<<s<<"\n";
		if (parameters.saveEach>0)
			os<<"parameters.SaveEach="<<parameters.saveEach<<"\n";
		os<<"parameters.latticeLength="<<parameters.latticeLength<<"\n";
		os<<"parameters.coresForKernel="<<parameters.coresForKernel<<"\n";
		os<<"parameters.adjustEach="<<parameters.adjustEach<<"\n";

		typedef typename PsimagLite::Map<PsimagLite::String,RealType>::Type MapType;
		typedef typename MapType::const_iterator MapIterator;
		const MapType& mymap = parameters.mcWindow;
		for (MapIterator iter = mymap.begin(); iter != mymap.end(); iter++)
			os<<"parameters.mcWindow["<<iter->first<<"]="<<iter->second<<"\n";

		return os;
	}
} // namespace Spf
/*@}*/

#endif
