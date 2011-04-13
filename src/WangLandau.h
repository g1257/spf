
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
		enum {THERMALIZATION,MEASUREMENT};

		void set(
				const RealType& minE,
				const RealType& maxE,
				size_t steps,
				const RealType& lf,
				size_t eachCountForF)
		{
			lg_.init(minE,maxE,steps,1.0);
			h_.init(minE,maxE,steps);
			lf_ = lf;
			eachCountForF_ = eachCountForF;
		}

		//! Do we accept or not
		bool operator()(const RealType& enew,const RealType& eold)
		{
			RealType p = random_.random();
			RealType m = exp(getG(eold))/exp(getG(enew));
			int chk = lg_.add(eold,lf_);
			if (chk==1) throw std::runtime_error("WangLandau: E out of range\n");
			h_.add(eold);
			if (p<m) return true;
			return false;
		}

		void changeF(size_t iter,size_t phase)
		{
			if (phase == THERMALIZATION) return;
			if (iter>0 && iter % eachCountForF_ == 0) {
				lf_ *= 0.5;
			}
		}

		const RealType& f() const { return lf_; }

		void print(std::ostream& os) const
		{
			os<<"#WangLandauGE\n";
			for (size_t i=0;i<lg_.size();i++)  {
				os<<lg_.coorX(i)<<" "<<lg_.coorY(i)<<"\n";
			}
			os<<"#WangLandauH\n";
			for (size_t i=0;i<h_.size();i++)  {
				os<<h_.coorX(i)<<" "<<h_.coorY(i)<<"\n";
			}
			os<<"#WangLandauf="<<lf_<<"\n";
		}

	private:

		RealType getG(const RealType& e) const
		{
			size_t i = lg_.getSubInterval(e);
			return lg_.coorY(i);
		}

		Random48Type random_;
		HistogramType lg_;
		HistogramType h_;
		RealType lf_;
		size_t eachCountForF_;
		
	}; // WangLandau


	//! FIXME: EMPTY FOR NOW
//	std::ostream& operator<<(std::ostream& os,const WangLandau& wl)
//	{
//		return os;
//	}
} // namespace Spf

/*@}*/
#endif // PHONONS_2ORB_H
