/*
Copyright (c) 2009-2013, UT-Battelle, LLC
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

/*! \file HubbardOneOrbitalFields.h
 *
 *  Wrapper around classical spin, to say that only spin is a MC variable
 *
 */
#ifndef HUBBARD_ONE_ORB_FIELDS_H
#define HUBBARD_ONE_ORB_FIELDS_H
#include "loki/Typelist.h"
#include "ContVarFiniteOperations.h"

namespace Spf {
	template<typename FieldType,typename GeometryType>
	class HubbardOneOrbitalFields {

	public:

		enum {CHARGE, MAG};

		typedef ContVarFiniteOperations<GeometryType,FieldType> ContVarFiniteOperationsType;
		typedef typename ContVarFiniteOperationsType::DynVarsType ContVarFiniteType;
		typedef typename ContVarFiniteType::PairRealType PairRealType;
		
		typedef LOKI_TYPELIST_2(ContVarFiniteOperationsType,ContVarFiniteOperationsType) OperationsList;

		template<typename SomeParamsType>
		HubbardOneOrbitalFields(size_t vol,const SomeParamsType& params)
		    : name_(2),
		      charge_(vol,params.dynvarsfile,0,PairRealType(0,2)),
		      mag_(vol,params.dynvarsfile,1,PairRealType(-1,1))
		{
			setNames();
		}
		
		HubbardOneOrbitalFields(const ContVarFiniteType& charge,const ContVarFiniteType& mag)
		    : name_(2),charge_(charge),mag_(mag)
		{
			setNames();
		}
		
		const PsimagLite::String& name(size_t i) const
		{
			return name_[i];
		}
		
		void getField(ContVarFiniteType const** field,size_t i) const
		{
			assert(i == 0 || i == 1);
			if (i == 0)
				*field = &charge_;
			else if (i == 1)
				*field = &mag_;
			else
				throw PsimagLite::RuntimeError("HubbardOneOrbitalFields::getField()\n");
		}

		void getField(ContVarFiniteType** field,size_t i)
		{
			assert(i == 0 || i == 1);
			if (i==0)
				*field = &charge_;
			else if (i == 1)
				*field = &mag_;
			else
				throw PsimagLite::RuntimeError("HubbardOneOrbitalFields::getField()\n");
		}

		const ContVarFiniteType& getField(size_t what) const
		{
			return (what==CHARGE) ? charge_ : mag_;
		}

		
		template<typename FieldType2,typename GeometryType2>
		friend std::ostream& operator<<(std::ostream& os,const HubbardOneOrbitalFields<FieldType2,GeometryType2>& f);
		
	private:

		void setNames()
		{
			assert(name_.size()==2);
			name_[CHARGE] = "Charge";
			name_[MAG] = "Mag";
		}

		PsimagLite::Vector<PsimagLite::String>::Type name_;
		ContVarFiniteType charge_;
		ContVarFiniteType mag_;
		
	}; // HubbardOneOrbitalFields
	
	template<typename FieldType,typename GeometryType>
	std::ostream& operator<<(std::ostream& os,const HubbardOneOrbitalFields<FieldType,GeometryType>& f)
	{
		os<<f.charge_;
		os<<f.mag_;
		return os;
	}

} // namespace Spf

/*@}*/
#endif
