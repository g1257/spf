/*
Copyright (c) 2009-2014, UT-Battelle, LLC
All rights reserved

[Spf, Version 7.0.0]
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

/*! \file ParametersDmsMultiOrbital.h
 *
 *   Parameters for the DmsMultiOrbital model
 *
 */
#ifndef PARAMS_DMS_MULTIORBITA_H
#define PARAMS_DMS_MULTIORBITA_H
#include <complex>
#include <iostream>
#include <vector>
#include <iostream>

namespace Spf {

// please don't add member functions, this is a struct!
template<typename ParametersEngineType,typename IoInType>
struct ParametersDmsMultiOrbital {

	typedef typename ParametersEngineType::RealType RealType;
	typedef std::complex<RealType> ComplexType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

	//! ctor to read Model Parameters from inp file.
	ParametersDmsMultiOrbital(IoInType& io,
	                          const ParametersEngineType& engineParams)
	{
		io.readline(orbitals,"Orbitals=");

		// legacy reading of pair of doubles to make a complex:
		typename PsimagLite::Vector<RealType>::Type tmpReal;
		io.read(tmpReal,"Hoppings");
		if (engineParams.options.find("realhoppings") != PsimagLite::String::npos) {
			hoppings.resize(tmpReal.size());
			for (SizeType i=0;i<hoppings.size();i++)
				hoppings[i] = tmpReal[i];
		} else {
			hoppings.resize(tmpReal.size()/2);
			for (SizeType i=0;i<hoppings.size();i++)
				hoppings[i] = ComplexType(tmpReal[2*i],tmpReal[2*i+1]);
		}

		io.readline(J,"CouplingJ=");

		io.read(potentialV,"PotentialV");
		SizeType n = potentialV.size();

		io.readline(jafNn,"PARAMETERSJ_AF=");

		io.readline(jafNnn,"PARAMETERSJ_AF_NN=");

		io.readline(magneticField,"MagneticField=");

		io.readline(spinOrbitCoupling,"SPIN_ORBIT_COUPLING=");

		try {
			io.read(dmNn,"PARAMETERS_DM_NN");
		} catch (std::exception& e) {}

		try {
			io.read(dmNnn,"PARAMETERS_DM_NNN");
		} catch (std::exception& e) {}

		PsimagLite::Vector<SizeType>::Type tmp;
		try {
			io.read(tmp,"MODULUS");
		} catch (std::exception& e) {
			tmp.resize(n);
			for (SizeType i=0;i<tmp.size();i++) tmp[i] = i;
		}

		modulus.resize(n);
		for (SizeType i=0;i<modulus.size();i++) modulus[i] = 0;
		for (SizeType i=0;i<tmp.size();i++) modulus[tmp[i]] = 1;

		io.read(histogramParams,"HISTOGRAM");
	}

	SizeType orbitals;

	// packed as gamma1+gamma2*dof + dir*4
	// where dir goes from 0 to 11
	// and dof = 2*orbitals
	// and gamma1 = orb + spin*ORBITALS I think
	// or is it spin + orb*2 ?
	typename PsimagLite::Vector<ComplexType>::Type hoppings;

	// J value
	RealType J;

	// Onsite potential values, one for each site
	typename PsimagLite::Vector<RealType>::Type potentialV;

	// JAF n-n
	RealType jafNn;

	// JAF n-n-n
	RealType jafNnn;

	// Magnetic Field B for Zeeman term
	RealType magneticField;

	RealType spinOrbitCoupling; // =0.34

	VectorRealType dmNn;
	VectorRealType dmNnn;

	// Modulus (FIXME: use less storage here it should be either 0 or 1)
	PsimagLite::Vector<SizeType>::Type modulus;

	typename PsimagLite::Vector<RealType>::Type histogramParams;
}; // struct ParametersDmsMultiOrbital

//! Function that prints model parameters to stream os
template<typename ParametersEngineType,typename IoInType>
std::ostream& operator<<(
        std::ostream &os,
        const ParametersDmsMultiOrbital<ParametersEngineType,IoInType>& parameters)
{
	//os<<"parameters.nOfElectrons="<<parameters.nOfElectrons<<"\n";
	os<<"parameters.jafNn="<<parameters.jafNn<<"\n";
	os<<"parameters.jafNnn="<<parameters.jafNnn<<"\n";
	os<<"parameters.J="<<parameters.J<<"\n";
	os<<"parameters.potentialV\n";
	os<<parameters.potentialV;
	os<<"parameters.hoppings\n";
	os<<parameters.hoppings;
	os<<"parameters.magneticField="<<parameters.magneticField<<"\n";
	os<<"parameters.spinOrbitCoupling="<<parameters.spinOrbitCoupling<<"\n";
	if (parameters.dmNn.size() > 0) os<<"parameters.dmNn="<<parameters.dmNn;
	if (parameters.dmNnn.size() > 0) os<<"parameters.dmNnn="<<parameters.dmNnn;

	os<<"modulus\n";
	for (SizeType i=0;i<parameters.modulus.size();i++)
		if (parameters.modulus[i]!=0) os<<i<<" ";
	os<<"\n";
	os<<"histogramParams\n";
	os<<parameters.histogramParams;
	return os;
}
} // namespace Spf

/*@}*/
#endif // PARAMS_DMS_MULTIORBITA_H

