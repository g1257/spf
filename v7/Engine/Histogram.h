
/** \ingroup SPF */
/*@{*/

/*! \file Histogram.h
 *
 *  Histogram to store I(omega) observables
 *
 */
#ifndef HISTOGRAM_H_
#define HISTOGRAM_H_
#include "Utils.h"
#include "ProgressIndicator.h"
#include "TypeToString.h"

namespace Spf {

	template<typename RealType,typename FieldTyp>
	class Histogram {
	public:	

		void Histogram(
				const RealType& minE,
				const RealType& maxE,
				size_t steps)
		: minE_(minE),maxE_(maxE),steps_(steps),
		  histX_(steps+1),histY_(steps+1)
		{
			checkBounds();

			RealType deltaE = (maxE-minE)/steps;

			for (size_t i=0;i<histE.size();i++) {
				histE[i]=minE_+i*deltaE;
				histDE[i]=0;
			}
		}

		void add(const RealType& x,const RealType& y)
		{
			int n;
			size_t n = size_t(steps_*(x-minE_));
			n = size_t(n/(maxE-minE));
			// Don't remove this checking it's very important!
			if (n>=steps || n<0) {
				std::string s = "Histogram::add(" + ttos(x) + "," +
						ttos(y) + ") out of range\n";
				throw std::runtime_error(s.c_str());
			}
			histY_[n] += y;
		}

		//! divide all energies by a constant factor
		void divide(const RealType& div1)
		{
			RealType div2=(maxE_-minE_)/steps_;
			divideInternal(div1*div2);
		}

	private:

		void divideInternal(const RealType& div1)
		{
			for (size_t i=0;i<histY.size();i++)
				histDE[i] /= div1;
		}


		void checkBounds() const
		{
			if (steps_>0 && minE_<maxE_) return;

			std::string s ="steps=" + ttos(steps) + " and minE=" +
					ttos(minE) + " and maxE=" + maxE +"\n";
			s += "Histogram: " + __FILE__ + ":" + __LINE__ + "\n";
			throw std::runtime_error(s.c_str());
		}

	
	}; // class Histogram

} // namespace Spf

/*@}*/
#endif // HISTOGRAM_H_
