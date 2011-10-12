
/** \ingroup SPF */
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

namespace Spf {

	class TpemSubspace {
		
		TpemSubspace(size_t size)
		: flags_(size,0)
		{
			//reset();
		}

// 		void reset()
// 		{
// 			t->top = t->stack;
// 		}

		size_t size() const { return stack_.size(); }

		size_t top() const { return stack_.top(); }

		void fill()
		{
			for (size_t state=0;state<flags_.size();state++) push(state);
		}

		void push (size_t state)
		{
			if (flags_[state] != 0) return;
			flags_[state] = 1;
			stack_.push(state);
		}
	private:
		
		std::vector<size_t> flags_;
		std::stack<size_t> stack_;
	}; // class TpemSupspace
} // namespace Spf

/*@}*/
#endif // TPEM_SUBSPACE_H
