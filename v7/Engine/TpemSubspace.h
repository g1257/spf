
/** \ingroup TPEM */
/*@{*/

/*! \file TpemSupspace.h
 * 
 *  
 *
 */
#ifndef TPEM_SUBSPACE_H
#define TPEM_SUBSPACE_H
#include <vector>
#include <stack>
#include <iostream>

namespace Tpem {

	template<typename RealType,typename SparseMatrixType>
	class TpemSubspace {
		typedef typename SparseMatrixType::value_type RealOrComplexType;
		
	public:
		TpemSubspace(size_t size,const RealType& eps)
		: eps_(eps),flags_(size,0)
		{
			//reset();
		}

		void clear()
		{
			stack_.clear();
		}

		size_t size() const { return stack_.size(); }

// 		size_t top() const { return stack_.top(); }

		const size_t& operator()(size_t x) const { return stack_[x]; }

		void fill()
		{
			for (size_t state=0;state<flags_.size();state++) push(state);
		}

		void push (size_t state)
		{
			if (flags_[state] != 0) return;
			flags_[state] = 1;
			stack_.push_back(state);
		}

		void sparseProduct(const SparseMatrixType& matrix,
		                   std::vector<RealOrComplexType>& dest,
		                   const std::vector<RealOrComplexType>& src)
		{
			for (size_t i = 0; i < matrix.rank(); i++) dest[i] = 0.0;
			
			size_t oldtop = size();
			
			/* loop over states that have been used so far */
			for (size_t p = 0; p <oldtop; p++) {
				size_t j = stack_[p];
				//if (real(src[j])==0 && imag(src[j])==0) continue;
				/* loop over nonzero elements of j^th row */
				for (int k = matrix.getRowPtr(j);k<matrix.getRowPtr(j+1);k++) {
					size_t i = matrix.getCol(k);
					RealOrComplexType t = src[j] * matrix.getValue(k);
					RealType u = std::real(t)*std::real(t) + std::imag(t)*std::imag(t);
					dest[i] += t;
					if (u > eps_) push(i);
				}
			}
		}
		
		static void sparseProductPem(const SparseMatrixType& matrix,
		                             std::vector<RealOrComplexType>& dest,
		                             const std::vector<RealOrComplexType>& src)
		{	
			for (size_t i = 0; i < matrix.rank(); i++) dest[i] = 0.0;

			/* loop over all rows */
			for (size_t j=0;j<matrix.rank();j++) {
				/* loop over nonzero elements of j^th row */	
				for (int k = matrix.getRowPtr(j);k<matrix.getRowPtr(j+1);k++) {
					size_t i = matrix.getCol(k);
					RealOrComplexType t = src[j] * matrix.getValue(k);
					dest[i] += t;
				}
			}
		}

	private:

		const RealType eps_;
		std::vector<size_t> flags_;
		std::vector<size_t> stack_;
	}; // class TpemSupspace
} // namespace Tpem

/*@}*/
#endif // TPEM_SUBSPACE_H
