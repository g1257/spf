// BEGIN LICENSE BLOCK
/*
Copyright © 2009 , UT-Battelle, LLC
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

/*! \file PhononsTwoOrbitalsFields.h
 *
 *  Wrapper around classical spin, to say that only spin is a MC variable
 *
 */
#ifndef PHONONS_2ORB_FIELDS_H
#define PHONONS_2ORB_FIELDS_H
#include "Spin.h"
#include "Phonon.h"

namespace Spf {
	template<typename FieldType>
	class PhononsTwoOrbitalsFields {
	public:
		typedef Spin<FieldType> SpinType;
		typedef Phonon<FieldType> PhononType;
		typedef SpinType Type0;
		typedef PhononType Type1;
		
		PhononsTwoOrbitalsFields(size_t vol,const std::string& mcstarttype) :
				spin_(vol,mcstarttype),phonon_(vol,mcstarttype)
		{}
		
		PhononsTwoOrbitalsFields(const SpinType& spin,const PhononType& phonon) : 
				spin_(spin),phonon_(phonon)
		{
		}
		
		size_t size() const { return 2; } //  spins and phonons for this model need MC simulation
		
		const std::string& name(size_t i) const { return name_[i]; }
		
		
// 		template<typename T>
// 		T& getField(T&)
// 		{
// 			throw std::runtime_error("Should not reach here!\n");
// 		}
		
		template<int num,typename SomeType>
		SomeType& getField();
		
		template<typename FieldType2>
		friend std::ostream& operator<<(std::ostream& os,PhononsTwoOrbitalsFields<FieldType2>& f);
		
	private:
		static std::vector<std::string> name_;
		SpinType spin_;
		PhononType phonon_;
		
	}; // PhononsTwoOrbitalsFields
	
	template<typename FieldType>
	std::ostream& operator<<(std::ostream& os,PhononsTwoOrbitalsFields<FieldType>& f)
	{
		os<<f.spin_;
		os<<f.phonon_;
		return os;
	}
	template<typename FieldType>
	std::vector<std::string> PhononsTwoOrbitalsFields<FieldType>::name_(2);
	
	template<>
	template<>
	PhononsTwoOrbitalsFields<double>::SpinType& PhononsTwoOrbitalsFields<double>::getField<0, PhononsTwoOrbitalsFields<double>::SpinType>()	
	{
		return spin_;
	}

	template<>
	template<>
	PhononsTwoOrbitalsFields<double>::PhononType& PhononsTwoOrbitalsFields<double>::getField<1,PhononsTwoOrbitalsFields<double>::PhononType>()	
	{
		return phonon_;
	}
} // namespace Spf

/*@}*/
#endif
