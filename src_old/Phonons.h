#ifndef PHONONS_H_
#define PHONONS_H_

#include <iostream>
#include <stdexcept>

template<typename ParametersType,typename GeometryType>
class Phonons {
public:
	
	Phonons(const ParametersType& parameters,const GeometryType& geometry) :
		ether_(parameters), geometry_(geometry)
		{}
		
// 	void init(const ParametersType* parameters,const GeometryType* geometry)
// 	{
// 		ether_ = parameters;
// 		geometry_ = geometry;
// 	}
	
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
	const ParametersType& ether_;
	const GeometryType& geometry_;
	
};


#endif
