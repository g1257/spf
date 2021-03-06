
/** \ingroup SPF */
/*@{*/

/*! \file ClassicalSpinOperations.h
 *
 *  Operations related to a vector of classical spins
 *
 */
#ifndef CLASSICAL_SPIN_OPS_H
#define CLASSICAL_SPIN_OPS_H
#include "Vector.h"
#include "Spin.h"
#include "../Engine/ProgramGlobals.h"

namespace Spf {

	template<typename GeometryType_,typename RealType>
	class ClassicalSpinOperations {

		typedef PsimagLite::Vector<RealType> VectorType;
		typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

		static const bool isingSpins_ = false; // FIXME: make it runtime option

	public:

		typedef Spin<RealType> SpinType;
		typedef GeometryType_ GeometryType;
		typedef SpinType DynVarsType;

		template<typename SomeParamsType>
		ClassicalSpinOperations(const GeometryType& geometry,
		                        const SomeParamsType& params)
		: geometry_(geometry),
		  dynVars2_(0,params)
		{
			ProgramGlobals::checkMcWindow(params.mcWindow,1);
			mcwindowPhi_ = params.mcWindow[1];
		    mcwindowTheta_ = params.mcWindow[0];
		}

		void set(DynVarsType& dynVars)
		{
			dynVars_=&dynVars;
		}

		//! How to sweep the lattice
		template<typename RandomNumberGeneratorType>
		SizeType proposeSite(SizeType i,RandomNumberGeneratorType&) const
		{
			return i; //<-- zig-zag horizontal
			// zig-zag vertical:
			/*SizeType l = geometry_.length();
			SizeType x = i % l;
			SizeType y = i / l;
			return y + x*l;*/
			// random:
			//return SizeType(rng()*geometry_.volume());

		}

		template<typename RandomNumberGeneratorType>
		void proposeChange(SizeType i,RandomNumberGeneratorType& rng)
		{
			RealType thetaOld = dynVars_->theta[i];
			RealType phiOld = dynVars_->phi[i];

			dynVars2_ = *dynVars_;

			propose_(thetaOld,phiOld,dynVars2_.theta[i],dynVars2_.phi[i],rng);
		}

		void makeChange(SizeType i,const RealType& eps1,const RealType& eps2)
		{
			dynVars2_ = *dynVars_;
			dynVars2_.theta[i] += eps1;
			dynVars2_.phi[i] += eps2;
		}

		const DynVarsType& dynVars2() const { return dynVars2_; }

		RealType deltaDirect(SizeType i,RealType coupling1,RealType coupling2) const
		{
			SizeType z = geometry_.z(1)/2;
			VectorRealType coupling1v(z,coupling1);
			return deltaDirect(i,coupling1v,coupling2);
		}

		RealType deltaDirect(SizeType i,
		                       const VectorRealType& coupling1v,
		                       RealType coupling2) const
		{
			RealType sum = dSDirect(*dynVars_,dynVars2_,i,coupling1v);
			sum += directExchange2(dynVars2_,coupling2)
										-directExchange2(*dynVars_,coupling2);
			return sum;
		}

		RealType deltaMagneticField(SizeType i, const RealType& B) const
		{
			RealType dx = cos(dynVars2_.theta[i]) - cos(dynVars_->theta[i]);
			return dx * B;
		}

		RealType deltaDmInteraction(SizeType,
		                            const VectorRealType& dmNn,
		                            const VectorRealType& dmNnn) const
		{
			return calcDmTerm(dynVars2_,dmNn,dmNnn) - calcDmTerm(*dynVars_,dmNn,dmNnn);
		}

		RealType sineUpdate(SizeType i) const
		{
			RealType sineupdate= sin(dynVars_->theta[i]);
			if (sineupdate!=0) {
				sineupdate = sin(dynVars2_.theta[i])/sineupdate;
			} else {
				sineupdate = 1.0;
			}
			return sineupdate;
		}

		void accept(SizeType i)
		{
			dynVars_->theta[i]=dynVars2_.theta[i];
			dynVars_->phi[i]=dynVars2_.phi[i];
		}

		RealType directExchange2(const DynVarsType& dynVars,RealType coupling) const
		{
			SizeType n = dynVars.theta.size();
			RealType dS = 0;

			for (SizeType i=0;i<n;i++) {
				RealType t1=dynVars.theta[i];
				RealType p1=dynVars.phi[i];
				RealType cost1 = cos(t1);
				RealType sint1 = sin(t1);
				RealType cosp1 = cos(p1);
				RealType sinp1 = sin(p1);
				for (SizeType k = 0; k<geometry_.z(2); k++){
					SizeType j=geometry_.neighbor(i,k,2).first; /**next nearest neighbor */
					RealType t2=dynVars.theta[j];
					RealType p2=dynVars.phi[j];
					RealType tmp = cost1*cos(t2)+sint1*sin(t2)*(cosp1*cos(p2)+sinp1*sin(p2));
					dS += tmp;
				}
			}

			return coupling*dS*0.5;
		}

		RealType calcDmTerm(const DynVarsType& dynVars,
		                    const VectorRealType& dmNn,
		                    const VectorRealType& dmNnn) const
		{
			return calcDmTerm(dynVars,dmNn,1) + calcDmTerm(dynVars,dmNnn,2);

		}

		RealType calcDmTerm(const DynVarsType& dynVars,
		                    const VectorRealType& D,
		                    SizeType neighbor) const
		{
			if (D.size() < 3) return 0;
			SizeType n = dynVars.theta.size();
			RealType sum = 0;

			for (SizeType i=0;i<n;i++) {
				RealType t1=dynVars.theta[i];
				RealType p1=dynVars.phi[i];
				RealType cost1 = cos(t1);
				RealType sint1 = sin(t1);
				RealType cosp1 = cos(p1);
				RealType sinp1 = sin(p1);
				for (SizeType k = 0; k<geometry_.z(neighbor); k++){
					SizeType j=geometry_.neighbor(i,k,neighbor).first;
					if (i >= j) continue;
					RealType t2=dynVars.theta[j];
					RealType p2=dynVars.phi[j];
					int sign = 1;
					if (neighbor == 2) {
						if (geometry_.isBoundary(i,j,1)) sign = -1;
					} else {
						if (isBoundary(i,j)) sign = -1;
					}

					sum += D[0]*sign*(sint1*sinp1*cos(t2)-sin(t2)*sin(p2)*cos(t1));
					sum += D[1]*sign*(sin(t2)*cos(p2)*cost1-sint1*cosp1*cos(t2));
					sum += D[2]*sign*sin(t2)*sin(p2)*sint1*cosp1;
					sum -= D[2]*sign*sin(t2)*cos(p2)*sint1*sinp1;
				}
			}

			return sum;
		}

		RealType calcSuperExchange(const DynVarsType& dynVars,
				                              const RealType& coupling)
					const
		{
			SizeType z = geometry_.z(1);
			VectorRealType coupling1v(z,coupling);
			return calcSuperExchange(dynVars,coupling1v);

		}

		RealType calcSuperExchange(const DynVarsType& dynVars,
		                           const VectorRealType& coupling)
			const
		{
			RealType sum = 0;
			for (SizeType i=0;i<geometry_.volume();i++) {
				for (SizeType k = 0; k<geometry_.z(1); k++){
					SizeType j=geometry_.neighbor(i,k).first;
					SizeType dir=geometry_.neighbor(i,k).second;

					RealType t1=dynVars.theta[i];
					RealType t2=dynVars.theta[j];
					RealType p1=dynVars.phi[i];
					RealType p2=dynVars.phi[j];
					RealType tmp = cos(t1)*cos(t2)+sin(t1)*sin(t2)*(cos(p1)*cos(p2)+sin(p1)*sin(p2));
					assert(dir<coupling.size());
					sum += coupling[dir]*tmp;
				}
			}
			return sum*0.5;
		}

		RealType calcMag(const DynVarsType& dynVars) const
		{
			VectorRealType mag(3);

			for (SizeType i=0;i<geometry_.volume();i++) {
				mag[0] += sin(dynVars.theta[i])*cos(dynVars.phi[i]);
				mag[1] += sin(dynVars.theta[i])*sin(dynVars.phi[i]);
				mag[2] += cos(dynVars.theta[i]);
			}
			return (mag[0]*mag[0]+mag[1]*mag[1]+mag[2]*mag[2]);
		}

		//! For diluted systems
		RealType calcMag(
				const DynVarsType& dynVars,
				const PsimagLite::Vector<SizeType>::Type& modulus) const
		{
			VectorRealType mag(3);

			for (SizeType i=0;i<geometry_.volume();i++) {
				if (modulus[i]==0) continue;
				mag[0] += sin(dynVars.theta[i])*cos(dynVars.phi[i]);
				mag[1] += sin(dynVars.theta[i])*sin(dynVars.phi[i]);
				mag[2] += cos(dynVars.theta[i]);
			}
			return (mag[0]*mag[0]+mag[1]*mag[1]+mag[2]*mag[2]);
		}

		void calcMagVector(
				VectorRealType& mag,
				const DynVarsType& dynVars) const
		{
			for (SizeType i=0;i<geometry_.volume();i++) {
				mag[0] += sin(dynVars.theta[i])*cos(dynVars.phi[i]);
				mag[1] += sin(dynVars.theta[i])*sin(dynVars.phi[i]);
				mag[2] += cos(dynVars.theta[i]);
			}
		}


		void classicalCorrelations(VectorType &cc,
				 //PsimagLite::Vector<RealType>::Type& weight,
				 const DynVarsType& dynVars)
		{
			SizeType n = geometry_.volume();

			for (SizeType i=0;i<n;i++) {
				RealType temp=0;
				SizeType counter=0;
				for (SizeType j=0;j<n;j++) {
					SizeType k = geometry_.add(i,j);
					//cerr<<"calcClasCor: "<<i<<"+"<<j<<"="<<k<<endl;
					//if (ether.modulus[k]==0 || ether.modulus[j]==0) continue;
					temp+= cos(dynVars.theta[k])*cos(dynVars.theta[j])+
							sin(dynVars.theta[k])*sin(dynVars.theta[j])*
							cos(dynVars.phi[k]-dynVars.phi[j]);
					counter++;
				}
				if (counter>0) temp /= counter;
				//weight[i]=counter;
				cc[i] += temp;
			}
		}


	private:

		bool isBoundary(SizeType i, SizeType j) const
		{
			SizeType dim = geometry_.dim();
			for (SizeType d = 0; d < dim; ++d) {
				if (geometry_.isBoundary(i,j,d)) return true;
			}

			return false;
		}

		template<typename RngType>
		void propose_(
				RealType thetaOld,
				RealType phiOld,
				RealType &thetaNew,
				RealType &phiNew,
				RngType& rng)
		{
			//if (fabs(mcwindowPhi_)<1e-8 && fabs(mcwindow_[1]<1e-8)) return;

			if (isingSpins_) {
				if (thetaOld==0) thetaNew=M_PI;
				else thetaNew=0;
				phiNew=0;
				return;
			}

			if (mcwindowTheta_ < 0) {
				thetaNew=2*rng()- 1;
			} else {
				thetaNew = thetaOld + mcwindowTheta_*(rng() - 0.5);
			}

			while (thetaNew < -1) thetaNew += 1;
			while (thetaNew > 1) thetaNew -= 1;
			assert(fabs(thetaNew)<1);
			thetaNew = acos(thetaNew);

			if (mcwindowPhi_<0) {
				phiNew = 2*M_PI*rng();
			} else {
				phiNew=phiOld+2*M_PI*(rng()- 0.5)*mcwindowPhi_;
			}
			/*if (ether.isSet("sineupdate")) {
				thetaNew = M_PI*rng();
			}*/

			while (thetaNew<0) {
				thetaNew = -thetaNew;
				phiNew+=M_PI;
			}
			while (thetaNew>M_PI) {
				thetaNew -= M_PI;
				phiNew+=M_PI;
			}

			while (phiNew<0) phiNew += 2*M_PI;
			while (phiNew>2*M_PI) phiNew -= 2*M_PI;
			//std::cerr<<"ThetaOld="<<thetaOld<<" thetaNew="<<thetaNew<<"\n";
			//std::cerr<<"PhiOld="<<phiOld<<" phiNew="<<phiNew<<"\n";
		}

		RealType dSDirect(const DynVarsType& dynVars,
		                    const DynVarsType& dynVars2,
		                    SizeType i,
		                    VectorRealType coupling) const
		{
			RealType dS = 0;

			for (SizeType k = 0; k<geometry_.z(1); k++){
				SizeType j=geometry_.neighbor(i,k).first;
				SizeType dir = geometry_.neighbor(i,k).second;
				assert(dir<coupling.size());
				RealType tmp = (sin(dynVars2.theta[i])*cos(dynVars2.phi[i])-sin(dynVars.theta[i])*
					cos(dynVars.phi[i]))*sin(dynVars.theta[j])*cos(dynVars.phi[j]) +
						(sin(dynVars2.theta[i])*sin(dynVars2.phi[i])-sin(dynVars.theta[i])*
					sin(dynVars.phi[i]))*sin(dynVars.theta[j])*sin(dynVars.phi[j]) +
					(cos(dynVars2.theta[i])-cos(dynVars.theta[i]))*cos(dynVars.theta[j]);
				dS += coupling[dir]*tmp;
			}
			//if (ether.isSet("magneticfield")) tmp = Zeeman(dynVars2,geometry,ether)-Zeeman(dynVars,geometry,ether);
			//else tmp =0;
			// dS += tmp;

			return dS;
		}

		const GeometryType& geometry_;
		RealType mcwindowPhi_;
		RealType mcwindowTheta_;
		DynVarsType* dynVars_;
		DynVarsType dynVars2_;


	}; // Engine
} // namespace Spf

/*@}*/
#endif // CLASSICAL_SPIN_OPS_H
