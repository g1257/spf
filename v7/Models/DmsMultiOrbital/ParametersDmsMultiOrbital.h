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

/*! \file ParametersDmsMultiOrbital.h
 *
 *   Parameters for the DmsMultiOrbital model
 *
 */
#ifndef PARAMS_DMS_MULTIORBITA_H
#define PARAMS_DMS_MULTIORBITA_H
#include "Utils.h"
#include "SimpleReader.h"

namespace Spf {
	template<typename Field>
	struct ParametersDmsMultiOrbital {
		// total number of sites in the system
		int linSize;
		

		// packed as gamma1+gamma2*dof + dir*4
		// where dir=0 is FIXME
		// and dof = 2*orbitals
		// and gamma1 = orb + spin*ORBITALS I think
		// or is it spin + orb*2 ?
		std::vector<Field> hoppings; 
		// J value
		Field J;
		// Onsite potential values, one for each site
		std::vector<Field> potentialV;
		
		// JAF n-n
		Field jafNn;
		
		// JAF n-n-n
		Field jafNnn;

		Field spinOrbitCoupling; // =0.34

		// Modulus (FIXME: use less storage here it should be either 0 or 1)
		std::vector<size_t> modulus;
	}; // struct ParametersDmsMultiOrbital

	//! Operator to read Model Parameters from inp file.
	template<typename FieldType>
	ParametersDmsMultiOrbital<FieldType>& operator<=(
			ParametersDmsMultiOrbital<FieldType>& parameters,
			Dmrg::SimpleReader& reader)
	{
		reader.read(parameters.linSize);
		reader.read(parameters.hoppings);
		
		reader.read(parameters.J);
		reader.read(parameters.potentialV);

		reader.read(parameters.jafNn);
		reader.read(parameters.jafNnn);
		
		reader.read(parameters.spinOrbitCoupling);

		std::vector<size_t> tmp;
		reader.read(tmp);
		parameters.modulus.resize(parameters.linSize);
		for (size_t i=0;i<parameters.modulus.size();i++)
			parameters.modulus[i] = 0;
		for (size_t i=0;i<tmp.size();i++) parameters.modulus[tmp[i]] = 1;

		return parameters;
	}
	
	//! Function that prints model parameters to stream os
	template<typename FieldType>
	std::ostream& operator<<(
			std::ostream &os,
			const ParametersDmsMultiOrbital<FieldType>& parameters)
	{
		os<<"parameters.linSize="<<parameters.linSize<<"\n";
		//os<<"parameters.nOfElectrons="<<parameters.nOfElectrons<<"\n";
		os<<"parameters.jafNn="<<parameters.jafNn<<"\n";
		os<<"parameters.jafNnn="<<parameters.jafNnn<<"\n";
		os<<"parameters.J="<<parameters.J<<"\n";
		os<<"parameters.potentialV\n";
		os<<parameters.potentialV;
		os<<"parameters.hoppings\n";
		os<<parameters.hoppings;
		os<<"parameters.spinOrbitCoupling="<<parameters.spinOrbitCoupling<<"\n";
		os<<"modulus\n";
		for (size_t i=0;i<parameters.modulus.size();i++)
			if (parameters.modulus[i]!=0) os<<i<<" ";
		os<<"\n";
		return os;
	}
} // namespace Spf

/*@}*/
#endif // PARAMS_DMS_MULTIORBITA_H
