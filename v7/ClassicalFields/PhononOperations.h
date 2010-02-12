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

namespace Spf {
	template<typename GeometryType,typename DynVarsType>
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
		
		template<typename RandomNumberGeneratorType>
		void propose(size_t i,RandomNumberGeneratorType& rng)
		{
			PhononType phononsOld = dynVars->phonons[i];
			
			dynVars2_ = *dynVars_;
			
			propose_(phononsOld,dynVars2_.phonons[i],rng);
		}
		
		const DynVarsType& dynVars2() const { return dynVars2_; } 
		
		FieldType deltaDirect(size_t i,FieldType coupling1,FieldType coupling2) const
		{
			FieldType sum = dSDirect(*dynVars_,dynVars2_,i,coupling1);
			sum += directExchange2(dynVars2_,coupling2)
						-directExchange2(*dynVars_,coupling2);
			return sum;
		}
		
				
		FieldType sineUpdate(size_t i) const
		{
			FieldType sineupdate= sin(dynVars_->theta[i]);
			if (sineupdate!=0) {
				sineupdate = sin(dynVars2_.theta[i])/sineupdate;
			} else {
				sineupdate = 1.0;
			}
			return sineupdate;
		}
		
		void accept(size_t i)
		{
			dynVars_->theta[i]=dynVars2_.theta[i];
			dynVars_->phi[i]=dynVars2_.phi[i];
		}
	
		
		template<typename DynVarsType>
		double calcPhononDiff(int direction,int ind,DynVarsType const &dynVars) const
		{
			if (direction >= geometry_.dim()) return 0; 
			int j = geometry_.neighbor(ind,2*direction+1);
			return  (dynVars.phonons[ind][direction]-dynVars.phonons[j][direction]);
		}
		
		template<typename DynVarsType>
		double calcPhonon(int ind,DynVarsType const &dynVars,int what) const
		{
			double ret=0;
			double sqrt3=1.732050807569;
			double sqrt2=1.414213562373;
			double sqrt6=2.449489742783;
			
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
				std::cerr<<"I don't know what to calculate at "<<__FILE__<<" "<<__LINE__<<" with what="<<what<<std::endl;
				throw std::runtime_error("Phonons class\n");
			}
			return ret;
		}
		
		template<typename FieldType,typename DynVarsType>
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
		
		template<typename RandomNumberGeneratorType>
		void propose_(const PhononType& phononsOld,PhononType& phononsNew,RandomNumberGeneratorType& rng)
		{
			for (size_t i=0;i<phononsNew.size();i++) {
				phononsNew[i]=phononsOld[i] + (rng()- 0.5)*mcwindow_;
				//if (fabs(phononsNew[i]) > ether.maxPhonons) phononsNew[i]= 0.9*ether.maxPhonons;
			}
		}

		FieldType dSDirect(const DynVarsType& dynVars,const DynVarsType& dynVars2, size_t i,FieldType coupling) const
		{
			double dS=0;

			for (size_t alpha=0;alpha<dynVars.phonons[i].size();alpha++) {
				tmp = square(calcPhonon(i,dynVars2,alpha))
				- square(calcPhonon(i,dynVars,alpha));
				for (size_t k=0;k<geometry.z(1);k++) {
					size_t j = geometry.neighbor(i,k).first;
					tmp += square(calcPhonon(j,dynVars2,alpha));
					tmp -= square(calcPhonon(j,dynVars,alpha));
				}
				dS += coupling[alpha]*tmp;
			}
			return dS;
		}
	};

} // namespace Spf

/*@}*/
#endif
