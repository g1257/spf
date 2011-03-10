/** \ingroup SPF */
/*@{*/

/*! \file ParametersPhononsTwoOrbitals.h
 *
 *  Contains the parameters for the ParametersPhononsTwoOrbitals  model and function to read them from a JSON file
 *
 */

#ifndef PHONON_OPS_H_
#define PHONON_OPS_H_

#include "Utils.h"
#include "Vector.h"

namespace Spf {
	template<typename GeometryType,typename DynVarsType> // DynVarsType == PhononsType
	class PhononOperations {
		typedef typename DynVarsType::FieldType FieldType;
		typedef std::vector<FieldType> PhononType;
		
	public:
		PhononOperations(const GeometryType& geometry,FieldType mcwindow) 
			: geometry_(geometry),mcwindow_(mcwindow),dynVars2_(0,"none")
		{
		}
		
		void set(DynVarsType& dynVars)
		{
			dynVars_=&dynVars;
		}

		//! How to sweep the lattice
		template<typename RngType>
		size_t proposeSite(size_t i,RngType& rng) const
		{
			return i; //<-- zig-zag horizontal
			// zig-zag vertical:
			/*size_t l = geometry_.length();
			size_t x = i % l;
			size_t y = i / l;
			return y + x*l;*/
			// random:
			//return size_t(rng()*geometry_.volume());
			
		}
		
		
		template<typename RngType>
		void proposeChange(size_t i,RngType& rng)
		{
			PhononType phononsOld = dynVars_->phonon[i];
			
			dynVars2_ = *dynVars_;
			
			propose_(phononsOld,dynVars2_.phonon[i],rng);
		}
		
		const DynVarsType& dynVars2() const { return dynVars2_; } 
		
		FieldType deltaDirect(size_t i,const PhononType& coupling) const
		{
			return dSDirect(*dynVars_,dynVars2_,i,coupling);
		}

		FieldType sineUpdate(size_t i) const
		{
			return 1.0; // no measure for phonons
		}
		
		void accept(size_t i)
		{
			dynVars_->phonons[i]=dynVars2_.phonons[i];
		}
	
		FieldType calcPhononDiff(size_t direction,size_t ind,const DynVarsType& dynVars) const
		{
			if (direction >= geometry_.dim()) return 0; 
			size_t j = geometry_.neighbor(ind,2*direction+1).first;
			return  (dynVars.phonon[ind][direction]-dynVars.phonon[j][direction]);
		}
		
		FieldType calcPhonon(size_t ind,const DynVarsType& dynVars,size_t what) const
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

		void calcQvector(std::vector<FieldType>& v,size_t p,const DynVarsType& dynVars) const
		{
			size_t numberOfNormalModes = 3;
			v.resize(numberOfNormalModes);
			for (size_t i=0;i<numberOfNormalModes;i++)
				v[i]=calcPhonon(p,dynVars,i);
		}
		
	private:
		const GeometryType& geometry_;
		const FieldType& mcwindow_;
		DynVarsType* dynVars_;
		DynVarsType dynVars2_;
		
		template<typename RngType>
		void propose_(
				const PhononType& phononsOld,
				PhononType& phononsNew,
				RngType& rng)
		{
			for (size_t i=0;i<phononsNew.size();i++) {
				phononsNew[i]=phononsOld[i] + (RngType::random()- 0.5)*mcwindow_;
				//if (fabs(phononsNew[i]) > ether.maxPhonons) phononsNew[i]= 0.9*ether.maxPhonons;
			}
		}

		FieldType dSDirect(const DynVarsType& dynVars,const DynVarsType& dynVars2, size_t i,
				  const PhononType& coupling) const
		{
			double dS=0;

			for (size_t alpha=0;alpha<dynVars.phonons[i].size();alpha++) {
				FieldType tmp = square(calcPhonon(i,dynVars2,alpha))
				- square(calcPhonon(i,dynVars,alpha));
				for (size_t k=0;k<geometry_.z(1);k++) {
					size_t j = geometry_.neighbor(i,k).first;
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
