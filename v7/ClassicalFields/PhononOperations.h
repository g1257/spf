/** \ingroup SPF */
/*@{*/

/*! \file ParametersPhononsTwoOrbitals.h
 *
 *  Contains the parameters for the ParametersPhononsTwoOrbitals  model and function to read them from a JSON file
 *
 */

#ifndef PHONON_OPS_H_
#define PHONON_OPS_H_

#include "Vector.h"
#include "Phonon.h"

namespace Spf {
	template<typename GeometryType,typename FieldType>
	class PhononOperations {
		
	public:

		typedef Phonon<FieldType> PhononType;
		typedef typename PhononType::OnePhononType OnePhononType;
		typedef PhononType DynVarsType;

		PhononOperations(const GeometryType& geometry,FieldType mcwindow) 
			: geometry_(geometry),mcwindow_(mcwindow),dynVars2_(0,"none")
		{
		}
		
		void set(PhononType& dynVars)
		{
			dynVars_=&dynVars;
		}

		//! How to sweep the lattice
		template<typename RngType>
		SizeType proposeSite(SizeType i,RngType& rng) const
		{
			return i; //<-- zig-zag horizontal
			// zig-zag vertical:
			/*SizeType l = geometry_.length();
			SizeType x = i % l;
			SizeType y = i / l;
			return y + x*l;*/
			// random:
			//return SizeType(rng()*geometry_.volume());
			
		}
		
		
		template<typename RngType>
		void proposeChange(SizeType i,RngType& rng)
		{
			OnePhononType phononsOld = dynVars_->phonon[i];
			
			dynVars2_ = *dynVars_;
			
			propose_(phononsOld,dynVars2_.phonon[i],rng);
		}
		
		const PhononType& dynVars2() const { return dynVars2_; }
		
		FieldType deltaDirect(SizeType i,const OnePhononType& coupling) const
		{
			return dSDirect(*dynVars_,dynVars2_,i,coupling);
		}

		FieldType sineUpdate(SizeType i) const
		{
			return 1.0; // no measure for phonons
		}
		
		void accept(SizeType i)
		{
			dynVars_->phonon[i]=dynVars2_.phonon[i];
		}
	
		FieldType calcPhononDiff(SizeType direction,SizeType ind,const PhononType& dynVars) const
		{
			if (direction >= geometry_.dim()) return 0; 
			SizeType j = geometry_.neighbor(ind,2*direction+1).first;
			return  (dynVars.phonon[ind][direction]-dynVars.phonon[j][direction]);
		}
		
		FieldType calcPhonon(SizeType ind,const PhononType& dynVars,SizeType what) const
		{
			FieldType ret=0;
			FieldType sqrt3=1.732050807569;
			FieldType sqrt2=1.414213562373;
			FieldType sqrt6=2.449489742783;
			
			if (what==0)  { /* calc q1 */
				ret = (calcPhononDiff(0,ind,dynVars) + calcPhononDiff(1,ind,dynVars) +
				calcPhononDiff(2,ind,dynVars));
				ret /= sqrt3;
			} else if (what==1) { /* calc q2 */
				ret = (calcPhononDiff(0,ind,dynVars)-calcPhononDiff(1,ind,dynVars));
				ret /= sqrt2;
			} else if (what==2) { /* calc q3 */
				ret = (2*calcPhononDiff(2,ind,dynVars)-calcPhononDiff(0,ind,dynVars)-calcPhononDiff(1,ind,dynVars));
				ret /= sqrt6;
			} else {
				throw std::runtime_error("Phonons class\n");
			}
			return ret;
		}

		void calcQvector(typename PsimagLite::Vector<FieldType>::Type& v,SizeType p,const PhononType& dynVars) const
		{
			SizeType numberOfNormalModes = 3;
			v.resize(numberOfNormalModes);
			for (SizeType i=0;i<numberOfNormalModes;i++)
				v[i]=calcPhonon(p,dynVars,i);
		}
		
	private:
		const GeometryType& geometry_;
		const FieldType& mcwindow_;
		PhononType* dynVars_;
		PhononType dynVars2_;
		
		template<typename RngType>
		void propose_(
				const OnePhononType& phononsOld,
				OnePhononType& phononsNew,
				RngType& rng)
		{
			for (SizeType i=0;i<phononsNew.size();i++) {
				phononsNew[i]=phononsOld[i] + (rng()- 0.5)*mcwindow_;
				//if (fabs(phononsNew[i]) > ether.maxPhonons) phononsNew[i]= 0.9*ether.maxPhonons;
			}
		}

		FieldType dSDirect(const PhononType& dynVars,const PhononType& dynVars2, SizeType i,
				  const OnePhononType& coupling) const
		{
			FieldType dS=0;

			for (SizeType alpha=0;alpha<dynVars.phonons[i].size();alpha++) {
				FieldType tmp = square(calcPhonon(i,dynVars2,alpha))
				- square(calcPhonon(i,dynVars,alpha));
				for (SizeType k=0;k<geometry_.z(1);k++) {
					SizeType j = geometry_.neighbor(i,k).first;
					tmp += square(calcPhonon(j,dynVars2,alpha));
					tmp -= square(calcPhonon(j,dynVars,alpha));
				}
				dS += coupling[alpha]*tmp;
			}
			return dS;
		}
	}; // PhononOperations

} // namespace Spf

/*@}*/
#endif
