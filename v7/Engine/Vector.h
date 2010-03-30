
/** \ingroup PsimagLite */
/*@{*/

/*! \file Vector.h
 *
 *  
 *
 */
#ifndef VECTOR_H_
#define VECTOR_H_

#include <iostream>
#include <vector>

namespace PsimagLite {
	template<typename FieldType>
	class Vector {
			typedef Vector<FieldType> ThisType;
		public:
			Vector(size_t n) : data_(n) {}
			Vector() { }
			
			ThisType& operator/=(const FieldType& d)
			{
				if (d==0) return *this;
				for (size_t i=0;i<data_.size();i++)
					data_[i] /= d;
				return *this;
			}

			void resize(size_t x) { data_.resize(x); }
			
			size_t size() const { return data_.size(); }
			
			FieldType operator[](size_t i) const
			{
				return data_[i];
			}
			
			FieldType& operator[](size_t i)
			{
				return data_[i];
			}
			
		private:
			std::vector<FieldType> data_;
	}; // Vector
	
	template<typename FieldType>
	std::ostream& operator<<(std::ostream& os,const Vector<FieldType>& v) 
	{
		os<<v.size()<<"\n";
		for (size_t i=0;i<v.size();i++) os<<v[i]<<"\n";
		return os;
	}
	
		
 	template<typename FieldType> // move somewhere else
 	std::istream& operator>>(std::istream& is,Vector<FieldType>& v)
 	{
		size_t ss=0;
		is>>ss;
		v.resize(ss);
 		for (size_t i=0;i<v.size();i++) is>>v[i];
 		return is;
 	}
};

/*@}*/
#endif // VECTOR_H_
