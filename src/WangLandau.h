
/** \ingroup SPF */
/*@{*/

/*! \file WangLandau.h
 *
 *  Implmementation of Wang LAndau for SPF models
 *
 */
#ifndef WANG_LANDAU_H
#define WANG_LANDAU_H
#include "Random48.h"
#include "histogram.h" // lowercase!!!

namespace Spf {

	template<typename RealType>
	class WangLandau {
		typedef PsimagLite::Random48<RealType> Random48Type;
		typedef Histogram HistogramType;

	public:
		void set(
				const RealType& minE,
				const RealType& maxE,
				size_t steps,
				const RealType& f)
		{
			g_.init(minE,maxE,steps,1.0);
			h_.init(minE,maxE,steps);
			f_ = f;
		}

		//! Do we accept or not
		bool operator()(const RealType& enew,const RealType& eold)
		{
			RealType p = random_.random();
			RealType m = getG(eold)/getG(enew);
			g_.multiply(eold,f_);
			h_.add(eold);
			if (p<m) return true;
			return false;
		}

	private:

		RealType getG(const RealType& e) const
		{
			size_t i = g_.getSubInterval(e);
			return g_.coorY(i);
		}

		Random48Type random_;
		HistogramType g_;
		HistogramType h_;
		RealType f_;
		
	}; // WangLandau


	//! FIXME: EMPTY FOR NOW
//	std::ostream& operator<<(std::ostream& os,const WangLandau& wl)
//	{
//		return os;
//	}
} // namespace Spf

/*@}*/
#endif // PHONONS_2ORB_H
