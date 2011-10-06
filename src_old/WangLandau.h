
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
#include "IoSimple.h"

namespace Spf {

	template<typename RealType>
	class WangLandau {
		typedef PsimagLite::Random48<RealType> Random48Type;
		typedef Histogram HistogramType;

	public:
		enum {THERMALIZATION,MEASUREMENT};

		WangLandau() : enabled_(false) { }

		void init(
				const RealType& minE,
				const RealType& maxE,
				size_t steps,
				const RealType& lf)
//				size_t eachCountForF)
		{
			lg_.init(minE,maxE,steps,1.0);
			h_.init(minE,maxE,steps);
			lf_ = lf;
//			eachCountForF_ = eachCountForF;
			enabled_ = true;
		}

		void init(const std::string& file,size_t level)
		{
			std::string tmpString="#WangLandauGE";
			lg_.read(file,tmpString,level);
			tmpString="#WangLandauH";
			h_.read(file,tmpString,level);
			tmpString="#WangLandauf=";
			PsimagLite::IoSimple::In io(file.c_str());
			io.readline(lf_,"#WangLandauf=");
			print(std::cerr);
			if (lf_<0) throw std::runtime_error("WL: While reading lf\n");
			enabled_ = true;
		}

		void print(std::ostream& of) const
		{
			std::string tmpString="#WangLandauGE";
			lg_.print(of,tmpString);
			tmpString="#WangLandauH";
			h_.print(of,tmpString);
			of<<"#WangLandauf="<<lf_<<"\n";
		}

		//! Do we accept or not
		bool operator()(const RealType& enew,const RealType& eold)
		{
			RealType p = random_.random();
			RealType m = exp(getG(eold))/exp(getG(enew));
			int chk = lg_.add(eold,lf_);
			if (chk==1) throw std::runtime_error("WangLandau: E out of range\n");
			// if (chk==1) return false; // look the other way!!
			h_.add(eold);
			if (p<m) return true;
			return false;
		}

//		void changeF(size_t iter,size_t phase)
//		{
//			if (!enabled_) return;
//			if (phase == THERMALIZATION) return;
//			if (iter>0 && iter % eachCountForF_ == 0) {
//				lf_ *= 0.5;
//			}
//		}

//		const RealType& f() const { return lf_; }

	private:

		RealType getG(const RealType& e) const
		{
			size_t i = lg_.getSubInterval(e);
			return lg_.coorY(i);
		}

		bool enabled_;
		Random48Type random_;
		HistogramType lg_;
		HistogramType h_;
		RealType lf_;
//		size_t eachCountForF_;
		
	}; // WangLandau


	//! FIXME: EMPTY FOR NOW
//	std::ostream& operator<<(std::ostream& os,const WangLandau& wl)
//	{
//		return os;
//	}
} // namespace Spf

/*@}*/
#endif // PHONONS_2ORB_H
