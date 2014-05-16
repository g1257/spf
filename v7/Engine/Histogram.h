
/** \ingroup SPF */
/*@{*/

/*! \file Histogram.h
 *
 *  Histogram to store I(omega) observables
 *
 */
#ifndef HISTOGRAM_H_
#define HISTOGRAM_H_
#include "ProgressIndicator.h"
#include "TypeToString.h"

namespace Spf {

	template<typename RealType,typename FieldType>
	class Histogram {
	public:	

		Histogram(const RealType& minE,
		          const RealType& maxE,
		          SizeType steps)
		: minE_(minE),maxE_(maxE),steps_(steps),
		  histX_(steps+1),histY_(steps+1)
		{
			checkBounds();

			RealType deltaE = (maxE-minE)/steps;

			for (SizeType i=0;i<histX_.size();i++) {
				histX_[i]=minE_+i*deltaE;
				histY_[i]=0.0;
			}
		}

		void add(const RealType& x,const FieldType& y)
		{
			SizeType n = SizeType(steps_*(x-minE_));
			n = SizeType(n/(maxE_-minE_));
			// Don't remove this checking it's very important!
			if (n>=steps_ || x<minE_) {
				PsimagLite::String s = "Histogram::add(" + ttos(x) + "," +
						ttos(y) + ") out of range\n";
				//throw PsimagLite::RuntimeError(s.c_str());
				std::cerr<<s;
				return;
			}
			histY_[n] += y;
		}

		//! divide all energies by a constant factor
		RealType xWidth() const
		{
			return (maxE_-minE_)/steps_;
		}

		void divide(const RealType& div1)
		{
			for (SizeType i=0;i<histY_.size();i++)
				histY_[i] /= div1;
		}
		
		void reduce()
		{
			PsimagLite::MPI::reduce(histY_);
		}

		const RealType& x(SizeType i) const  { return histX_[i]; }

		const FieldType& y(SizeType i) const { return histY_[i]; }

		SizeType size() const { return steps_; }

	private:

		void checkBounds() const
		{
			if (steps_>0 && minE_<maxE_) return;

			PsimagLite::String s ="steps=" + ttos(steps_) + " and minE=" +
					ttos(minE_) + " and maxE=" + ttos(maxE_)+"\n";
			s += "Histogram: " + PsimagLite::String(__FILE__) + ":" +
					ttos(__LINE__) + "\n";
			throw PsimagLite::RuntimeError(s.c_str());
		}

		RealType minE_,maxE_;
		SizeType steps_;
		typename PsimagLite::Vector<RealType>::Type histX_;
		typename PsimagLite::Vector<FieldType>::Type histY_;
	
	}; // class Histogram

} // namespace Spf

/*@}*/
#endif // HISTOGRAM_H_
