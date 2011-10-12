
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
	public:
		TpemSubspace(size_t size)
		: flags_(size,0)
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
	private:
		
		std::vector<size_t> flags_;
		std::vector<size_t> stack_;
	}; // class TpemSupspace
} // namespace Spf

/*@}*/
#endif // TPEM_SUBSPACE_H
