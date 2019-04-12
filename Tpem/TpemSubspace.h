
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
		TpemSubspace(size_t size) : flags_(size,0),stack_(size,0),top_(0)
		{}

		void clear()
		{
			top_=0;
			for (size_t i=0;i<flags_.size();i++) flags_[i] = 0;
		}

// 		size_t size() const { return stack_.size(); }

 		const size_t& top() const { return top_; }

		const size_t& operator()(size_t x) const { return stack_[x]; }

		void fill()
		{
			for (size_t state=0;state<flags_.size();state++) push(state);
		}

		void push (size_t state)
		{
			if (flags_[state]) return;
			flags_[state] = 1;
			stack_[top_] = state;
			top_++;
		}

		void sparseProduct(const SparseMatrixType& matrix,
		                   std::vector<RealOrComplexType>& dest,
		                   const std::vector<RealOrComplexType>& src,
		                   const RealType& eps)
		{
			assert(matrix.rows() == matrix.cols());
			for (size_t i = 0; i < matrix.rows(); i++) dest[i] = 0.0;
			
			size_t oldtop = top_;
			
			/* loop over states that have been used so far */
			for (size_t p = 0; p <oldtop; p++) {
				size_t j = stack_[p];
				//if (real(src[j])==0 && imag(src[j])==0) continue;
				/* loop over nonzero elements of j^th row */
				int start = matrix.getRowPtr(j);
				int end = matrix.getRowPtr(j+1);
				for (int k = start;k<end;k++) {
					size_t i = matrix.getCol(k);
					RealOrComplexType t = src[j] * matrix.getValue(k);
					RealType u = std::real(t)*std::real(t) + std::imag(t)*std::imag(t);
					dest[i] += t;
					if (u > eps) push(i);
				}
			}
		}
		
		static void sparseProductPem(const SparseMatrixType& matrix,
		                             std::vector<RealOrComplexType>& dest,
		                             const std::vector<RealOrComplexType>& src)
		{	
			assert(matrix.rows() == matrix.cols());
			for (size_t i = 0; i < matrix.cols(); i++) dest[i] = 0.0;

			/* loop over all rows */
			for (size_t j=0;j<matrix.rows();j++) {
				/* loop over nonzero elements of j^th row */	
				for (int k = matrix.getRowPtr(j);k<matrix.getRowPtr(j+1);k++) {
					size_t i = matrix.getCol(k);
					RealOrComplexType t = src[j] * matrix.getValue(k);
					dest[i] += t;
				}
			}
		}

	private:
		std::vector<size_t> flags_;
		std::vector<size_t> stack_;
		size_t top_;
	}; // class TpemSupspace
} // namespace Tpem

/*@}*/
#endif // TPEM_SUBSPACE_H
