
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

namespace Spf {
	
	template<typename GeometryType_,typename DynVarsType_>
	class ClassicalSpinOperations {
		typedef typename DynVarsType_::FieldType RealType;
		typedef PsimagLite::Vector<RealType> VectorType;
			
		static const bool isingSpins_ = false; // FIXME: make it runtime option
		
	public:
		typedef GeometryType_ GeometryType;
		typedef DynVarsType_ DynVarsType;
		
		template<typename SomeParamsType>
		ClassicalSpinOperations(const GeometryType& geometry,const SomeParamsType& params) 
		: geometry_(geometry),mcwindow_(params.mcWindow),dynVars2_(0,params)
		{
		}
		
		void set(DynVarsType& dynVars)
		{
			dynVars_=&dynVars;
		}

		//! How to sweep the lattice
		template<typename RandomNumberGeneratorType>
		size_t proposeSite(size_t i,RandomNumberGeneratorType& rng) const
		{
			return i; //<-- zig-zag horizontal
			// zig-zag vertical:
			/*size_t l = geometry_.length();
			size_t x = i % l;
			size_t y = i / l;
			return y + x*l;*/
			// random:
			//return size_t(rng()*geometry_.volume());
			
		}
		
		template<typename RandomNumberGeneratorType>
		void proposeChange(size_t i,RandomNumberGeneratorType& rng)
		{
			RealType thetaOld = dynVars_->theta[i];
			RealType phiOld = dynVars_->phi[i];
			
			dynVars2_ = *dynVars_;
			
			propose_(thetaOld,phiOld,dynVars2_.theta[i],dynVars2_.phi[i],rng);
		}
		
		void makeChange(size_t i,const RealType& eps1,const RealType& eps2)
		{
			dynVars2_ = *dynVars_;
			dynVars2_.theta[i] += eps1;
			dynVars2_.phi[i] += eps2;
		}
		
		const DynVarsType& dynVars2() const { return dynVars2_; } 
		
		RealType deltaDirect(size_t i,RealType coupling1,RealType coupling2) const
		{
			size_t z = geometry_.z(1)/2;
			std::vector<RealType> coupling1v(z,coupling1);
			return deltaDirect(i,coupling1v,coupling2);
		}
		
		RealType deltaDirect(size_t i,
		                       const std::vector<RealType>& coupling1v,
		                       RealType coupling2) const
		{
			RealType sum = dSDirect(*dynVars_,dynVars2_,i,coupling1v);
			sum += directExchange2(dynVars2_,coupling2)
										-directExchange2(*dynVars_,coupling2);
			return sum;
		}

		RealType deltaMagneticField(size_t i, const RealType& B) const
		{
			RealType dx = cos(dynVars2_.theta[i]) - cos(dynVars_->theta[i]);
			return dx * B;
		}
				
		RealType sineUpdate(size_t i) const
		{
			RealType sineupdate= sin(dynVars_->theta[i]);
			if (sineupdate!=0) {
				sineupdate = sin(dynVars2_.theta[i])/sineupdate;
			} else {
				sineupdate = 1.0;
			}
			return sineupdate;
		}
		
		void accept(size_t i)
		{
			dynVars_->theta[i]=dynVars2_.theta[i];
			dynVars_->phi[i]=dynVars2_.phi[i];
		}
		
		RealType directExchange2(const DynVarsType& dynVars,RealType coupling) const
		{
			size_t n = dynVars.theta.size();
			RealType dS = 0;
			
			for (size_t i=0;i<n;i++) {
				RealType t1=dynVars.theta[i];
				RealType p1=dynVars.phi[i];
				RealType cost1 = cos(t1);
				RealType sint1 = sin(t1);
				RealType cosp1 = cos(p1);
				RealType sinp1 = sin(p1);
				for (size_t k = 0; k<geometry_.z(2); k++){
					size_t j=geometry_.neighbor(i,k,2).first; /**next nearest neighbor */
					RealType t2=dynVars.theta[j];
					RealType p2=dynVars.phi[j];
					RealType tmp = cost1*cos(t2)+sint1*sin(t2)*(cosp1*cos(p2)+sinp1*sin(p2));
					dS += tmp; 
				}
			}
			
			return coupling*dS*0.5;
		}
		RealType calcSuperExchange(const DynVarsType& dynVars,
				                              const RealType& coupling)
					const
		{
			size_t z = geometry_.z(1)/2;
			std::vector<RealType> coupling1v(z,coupling);
			return calcSuperExchange(dynVars,coupling1v);

		}
		
		RealType calcSuperExchange(const DynVarsType& dynVars,
		                              const std::vector<RealType>& coupling)
			const
		{
			RealType sum = 0;
			for (size_t i=0;i<geometry_.volume();i++) {
				for (size_t k = 0; k<geometry_.z(1); k++){
					size_t j=geometry_.neighbor(i,k).first;
					size_t dir=geometry_.neighbor(i,k).second;

					RealType t1=dynVars.theta[i];
					RealType t2=dynVars.theta[j];
					RealType p1=dynVars.phi[i];
					RealType p2=dynVars.phi[j];
					RealType tmp = cos(t1)*cos(t2)+sin(t1)*sin(t2)*(cos(p1)*cos(p2)+sin(p1)*sin(p2));
					sum += coupling[dir]*tmp;
				}
			}
			return sum*0.5;
		}

		RealType calcMag(const DynVarsType& dynVars) const
		{
			std::vector<RealType> mag(3);
			
			for (size_t i=0;i<geometry_.volume();i++) {
				mag[0] += sin(dynVars.theta[i])*cos(dynVars.phi[i]);
				mag[1] += sin(dynVars.theta[i])*sin(dynVars.phi[i]);
				mag[2] += cos(dynVars.theta[i]);
			}
			return (mag[0]*mag[0]+mag[1]*mag[1]+mag[2]*mag[2]);
		}

		//! For diluted systems
		RealType calcMag(
				const DynVarsType& dynVars,
				const std::vector<size_t>& modulus) const
		{
			std::vector<RealType> mag(3);

			for (size_t i=0;i<geometry_.volume();i++) {
				if (modulus[i]==0) continue;
				mag[0] += sin(dynVars.theta[i])*cos(dynVars.phi[i]);
				mag[1] += sin(dynVars.theta[i])*sin(dynVars.phi[i]);
				mag[2] += cos(dynVars.theta[i]);
			}
			return (mag[0]*mag[0]+mag[1]*mag[1]+mag[2]*mag[2]);
		}

		void calcMagVector(
				std::vector<RealType>& mag,
				const DynVarsType& dynVars) const
		{
			for (size_t i=0;i<geometry_.volume();i++) {
				mag[0] += sin(dynVars.theta[i])*cos(dynVars.phi[i]);
				mag[1] += sin(dynVars.theta[i])*sin(dynVars.phi[i]);
				mag[2] += cos(dynVars.theta[i]);
			}
		}

		
		void classicalCorrelations(VectorType &cc,
				 //std::vector<RealType> &weight,
				 const DynVarsType& dynVars)
		{
			size_t n = geometry_.volume();
			
			for (size_t i=0;i<n;i++) {
				RealType temp=0;
				size_t counter=0;
				for (size_t j=0;j<n;j++) {
					size_t k = geometry_.add(i,j);
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
		
		template<typename RngType>
		void propose_(
				RealType thetaOld,
				RealType phiOld,
				RealType &thetaNew,
				RealType &phiNew,
				RngType& rng)
		{
			if (fabs(mcwindow_[0])<1e-8 && fabs(mcwindow_[1]<1e-8)) return;

			if (isingSpins_) {
				if (thetaOld==0) thetaNew=M_PI; 
				else thetaNew=0;
				phiNew=0;
				return;
			} 
		
			if (mcwindow_[0]<0) {
				thetaNew = 2*rng.random()-1;
				phiNew = 2*M_PI*rng.random();
				thetaNew = acos(thetaNew);
			} else {
				thetaNew=2*rng.random()- 1;
				if (thetaNew < -1) thetaNew= 0;
				if (thetaNew > 1) thetaNew = 0;		
				phiNew=phiOld+2*M_PI*(rng.random()- 0.5)*mcwindow_[1];
				thetaNew = acos(thetaNew);
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
		                    size_t i,
		                    std::vector<RealType> coupling) const
		{
			RealType dS = 0;
				
			for (size_t k = 0; k<geometry_.z(1); k++){
				size_t j=geometry_.neighbor(i,k).first;
				size_t dir = geometry_.neighbor(i,k).second;
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
		const std::vector<RealType>& mcwindow_;
		DynVarsType* dynVars_;
		DynVarsType dynVars2_;
		
		
	}; // Engine
} // namespace Spf

/*@}*/
#endif // CLASSICAL_SPIN_OPS_H
