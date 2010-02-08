
/** \ingroup SPF */
/*@{*/

/*! \file PnictidesTwoOrbitals.h
 *
 *  PnictidesTwoOrbitals model
 *
 */
#ifndef PNICTIDES_2ORB_H
#define PNICTIDES_2ORB_H
#include "Utils.h"
#include "DynVars.h"
#include "RandomNumberGenerator.h"
#include "ProgressIndicator.h"

namespace Spf {
	template<typename FieldType,typename EngineParamsType,typename ParametersModelType,typename GeometryType>
	class PnictidesTwoOrbitals {
		
		typedef std::complex<FieldType> ComplexType;
		typedef psimag::Matrix<ComplexType> MatrixType;
		typedef RandomNumberGenerator<FieldType> RandomNumberGeneratorType;
		typedef typename GeometryType::PairType PairType;
		typedef Dmrg::ProgressIndicator ProgressIndicatorType;
		
		static const bool isingSpins_ = false;
		static const size_t nbands_ = 2;
		
		public:
			
		typedef DynVars<FieldType> DynVarsType;
		
		PnictidesTwoOrbitals(const EngineParamsType& engineParams,const ParametersModelType& mp,const GeometryType& geometry) :
			engineParams_(engineParams),mp_(mp),geometry_(geometry),hilbertSize_(2*nbands_*geometry_.volume()),
				      matrix_(hilbertSize_,hilbertSize_),progress_("PnictidesTwoOrbitals",0)
		{
		}
		
		size_t totalFlips() const { return geometry_.volume(); }
			
		size_t doMonteCarlo(DynVarsType& dynVars, size_t iter)
		{
			size_t acc = 0;
			std::vector<FieldType> eigOld;
			fillAndDiag(eigOld,dynVars);
			size_t n = geometry_.volume();
			
			for (size_t i=0;i<n;i++) {
				DynVarsType dynVars2 = dynVars;
				FieldType thetaNew,phiNew;
					
				r_newSpins(dynVars.theta[i],dynVars.phi[i],thetaNew,phiNew);
				
				dynVars2.theta[i]=thetaNew;
				dynVars2.phi[i]=phiNew;
				
				FieldType dsDirect = dSDirect(dynVars,dynVars2,i);
				dsDirect += directExchange2(dynVars2)-directExchange2(dynVars);
				
				//FieldType oldmu=engineParams_.mu;
				std::vector<FieldType> eigNew;
				fillAndDiag(eigNew,dynVars2);
				// if (ether.carriers>0) adjChemPot(eigNewAllBands); FIXME
				FieldType sineupdate= sin(dynVars.theta[i]);
				if (sineupdate!=0) {
					sineupdate = sin(dynVars2.theta[i])/sineupdate;
				} else {
					sineupdate = 1.0;
				}
				bool flag=doMetropolis(eigNew,eigOld,dsDirect,sineupdate);
					
				if (flag && !dynVars.isFrozen) { // Accepted
					dynVars.theta[i]=thetaNew;
					dynVars.phi[i]=phiNew;
					eigOld = eigNew;
					acc++;
				} else { // not accepted
					//engineParams_.mu=oldmu;
				}
			} // lattice sweep
			return acc;
		}
		
		void doMeasurements(DynVarsType& dynVars, size_t iter,std::ostream& fout)
		{
			std::string s = "iter=" + utils::ttos(iter); 
			
			progress_.printline(s,fout);
				
			std::vector<FieldType> eigs;
			fillAndDiag(eigs,dynVars);
				
			FieldType temp=calcNumber(eigs);
				
			s ="Number_Of_Electrons="+utils::ttos(temp);
			progress_.printline(s,fout);
			
			//s = "rankGlobal=";
			
			temp=calcElectronicEnergy(eigs);
			
			s="Electronic Energy="+utils::ttos(temp);
			progress_.printline(s,fout);
			
			FieldType temp2=calcSuperExchange(dynVars);
			s="Superexchange="+utils::ttos(temp2);
			progress_.printline(s,fout);
			
			temp += temp2;
			if (mp_.jafNnn!=0) {
				temp2=directExchange2(dynVars);
				s="Superexchange2="+utils::ttos(temp2);
				progress_.printline(s,fout);
				temp += temp2;
			}

			// total energy = electronic energy + superexchange + phonon energy
			s="TotalEnergy="+utils::ttos(temp);
			progress_.printline(s,fout);
				
			//s="Action=";
			
			//s="Number_Of_Holes=";
			
			//s="ChemPot="+utils::ttos(mu);
			//progress_.printline(s,fout);
			
			temp = calcMag(dynVars);
			s="Mag2="+utils::ttos(temp);
			progress_.printline(s,fout);
			
			temp=calcKinetic(dynVars,eigs);
			s ="KineticEnergy="+utils::ttos(temp);
			progress_.printline(s,fout);
		} // doMeasurements

		template<typename FieldType2,typename EngineParamsType2,typename ParametersModelType2,typename GeometryType2>
		friend std::ostream& operator<<(std::ostream& os,
				const PnictidesTwoOrbitals<FieldType2,EngineParamsType2,ParametersModelType2,GeometryType2>& model);
		
		private:
			
		void r_newSpins(FieldType thetaOld, FieldType phiOld, FieldType &thetaNew,FieldType &phiNew)
		{
			if (isingSpins_) {
				if (thetaOld==0) thetaNew=M_PI; 
				else thetaNew=0;
				phiNew=0;
			} else {
				if (engineParams_.mcWindow<0) {
					thetaNew = 2*rng_()-1;
					phiNew = 2*M_PI*rng_();
					thetaNew = acos(thetaNew);
				} else {
					thetaNew=2*rng_()- 1;
					if (thetaNew < -1) thetaNew= 0;
					if (thetaNew > 1) thetaNew = 0;		
					phiNew=phiOld+2*M_PI*(rng_()- 0.5)*engineParams_.mcWindow;
					thetaNew = acos(thetaNew);
				}
				/*if (ether.isSet("sineupdate")) {
					thetaNew = M_PI*ether.rng.myRandom();
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
			}
		}
		
		FieldType directExchange2(const DynVarsType& dynVars)
		{
			size_t n = geometry_.volume();
			FieldType dS = 0;
			
			for (size_t i=0;i<n;i++) {
				for (size_t k = 0; k<geometry_.z(2); k++){
					size_t j=geometry_.neighbor(i,k,2).first; /**next nearest neighbor */
					FieldType t1=dynVars.theta[i];
					FieldType t2=dynVars.theta[j];
					FieldType p1=dynVars.phi[i];
					FieldType p2=dynVars.phi[j];
					FieldType tmp = cos(t1)*cos(t2)+sin(t1)*sin(t2)*(cos(p1)*cos(p2)+sin(p1)*sin(p2));
					dS += mp_.jafNnn*tmp; 
				}
			}
			
			return dS*0.5;
		}
		
		FieldType dSDirect(const DynVarsType& dynVars,const DynVarsType& dynVars2, size_t i)
		{
			FieldType dS = 0;
				
			for (size_t k = 0; k<geometry_.z(1); k++){
				size_t j=geometry_.neighbor(i,k).first;
				FieldType tmp = (sin(dynVars2.theta[i])*cos(dynVars2.phi[i])-sin(dynVars.theta[i])*
					cos(dynVars.phi[i]))*sin(dynVars.theta[j])*cos(dynVars.phi[j]) +
						(sin(dynVars2.theta[i])*sin(dynVars2.phi[i])-sin(dynVars.theta[i])*
					sin(dynVars.phi[i]))*sin(dynVars.theta[j])*sin(dynVars.phi[j]) +
					(cos(dynVars2.theta[i])-cos(dynVars.theta[i]))*cos(dynVars.theta[j]);
				dS += mp_.jafNn*tmp;
			}
			//if (ether.isSet("magneticfield")) tmp = Zeeman(dynVars2,geometry,ether)-Zeeman(dynVars,geometry,ether);
			//else tmp =0;
			// dS += tmp;
		
			return dS;
		}
		
		void fillAndDiag(std::vector<FieldType> &eig,const DynVarsType& dynVars,char jobz='N')
		{
			if (jobz=='v') jobz='V';
			
			createHamiltonian(dynVars,matrix_);
			utils::diag(matrix_,eig,jobz);

			if (jobz!='V') sort(eig.begin(), eig.end(), std::less<FieldType>());
		}
		
		bool doMetropolis(const std::vector<FieldType>& eNew,const std::vector<FieldType>& eOld,
			FieldType dsDirect,FieldType sineupdate)
		{
			FieldType mu=engineParams_.mu;
			FieldType beta = engineParams_.beta;
			FieldType X =1.0;
			
			for (size_t i=0;i<eNew.size();i++) {
				FieldType temp = 0;
				if (eNew[i]>mu)
					temp = (double)(1.0+exp(-beta*(eNew[i]-mu)))/(1.0+exp(-beta*(eOld[i]-mu)));
				else
				temp =(double)(1.0+exp(beta*(eNew[i]-mu)))/
							(exp(beta*(eNew[i]-mu))+exp(-beta*(eOld[i]-eNew[i])));
			
				X *= temp;
			}
			
			//if (ether.isSet("sineupdate")) X *= sineupdate;
			X *=  exp(-beta*dsDirect);
			X = X/(1.0+X);
			
			FieldType r=rng_();
			
			if (X>r) return true;
			else return false;
		}
		
		void createHamiltonian (const DynVarsType& dynVars,MatrixType& matrix)
		{
			size_t volume = geometry_.volume();
			size_t norb = nbands_;
			size_t dof = norb * 2; // the 2 comes because of the spin
			std::vector<ComplexType> jmatrix(dof*dof);
			
			for (size_t gamma1=0;gamma1<matrix.n_row();gamma1++) 
				for (size_t p = 0; p < matrix.n_col(); p++) 
					matrix(gamma1,p)=0;
			
			for (size_t gamma1=0;gamma1<dof;gamma1++) {
				for (size_t p = 0; p < volume; p++) {
					
					auxCreateJmatrix(jmatrix,dynVars,p);
		
					size_t spin1 = size_t(gamma1/2);
					size_t orb1 = gamma1 % 2;
					matrix(p+gamma1*volume,p+gamma1*volume) = real(jmatrix[spin1+2*spin1]) + mp_.potentialV[p];
						
					for (size_t j = 0; j <  geometry_.z(1); j++) {	
						if (j%2!=0) continue;	
						PairType tmpPair = geometry_.neighbor(p,j);
						size_t k = tmpPair.first;
						size_t dir = tmpPair.second; // int(j/2);
						for (size_t orb2=0;orb2<norb;orb2++) {
							size_t gamma2 = orb2+spin1*norb; // spin2 == spin1 here
							matrix(p+gamma1*volume,k+gamma2*volume)=
								mp_.hoppings[orb1+orb2*norb+norb*norb*dir];
							matrix(k+gamma2*volume,p+gamma1*volume) = conj(matrix(p+gamma1*volume,k+gamma2*volume));
						}
					}
					//if (geometry.z(p,2)!=4 || geometry.z(p)!=4) throw std::runtime_error("neighbours wrong\n");
					for (size_t j = 0; j <  geometry_.z(2); j++) {
						if (j%2!=0) continue;	
						PairType tmpPair = geometry_.neighbor(p,j,2);
						size_t k = tmpPair.first;
						size_t dir = tmpPair.second;
						//std::cerr<<"Neigbors "<<p<<" "<<k<<"\n";
						for (size_t orb2=0;orb2<norb;orb2++) {
							size_t gamma2 = orb2+spin1*norb; // spin2 == spin1 here
							matrix(p+gamma1*volume,k+gamma2*volume)=
								mp_.hoppings[orb1+orb2*norb+norb*norb*dir];
							matrix(k+gamma2*volume,p+gamma1*volume) = conj(matrix(p+gamma1*volume,k+gamma2*volume));
						}
					}
					
					for (size_t spin2=0;spin2<2;spin2++) {
						if (spin1==spin2) continue; // diagonal term already taken into account
						size_t gamma2 = orb1+spin2*norb; // orb2 == orb1 here
						matrix(p+gamma1*volume,p + gamma2*volume)=jmatrix[spin1+2*spin2];
						matrix(p + gamma2*volume,p+gamma1*volume) = conj(matrix(p+gamma1*volume,p + gamma2*volume));
					}
				}
			}
		}
		
		void auxCreateJmatrix(std::vector<ComplexType>& jmatrix,const DynVarsType& dynVars,size_t site)
		{
			
			jmatrix[0]=0.5*cos(dynVars.theta[site]);
		
			jmatrix[1]=ComplexType(sin(dynVars.theta[site])*cos(dynVars.phi[site]),
				sin(dynVars.theta[site])*sin(dynVars.phi[site]));
		
			jmatrix[2]=conj(jmatrix[1]);
		
			jmatrix[3]= -cos(dynVars.theta[site]);
		
			for (size_t i=0;i<jmatrix.size();i++) jmatrix[i] *= mp_.J; 
		}
		
		FieldType calcNumber(const std::vector<FieldType>& eigs) const
		{
			FieldType sum=0;
			for (size_t i=0;i<eigs.size();i++) {
				sum += utils::fermi((eigs[i]-engineParams_.mu)*engineParams_.beta);
			}
			return sum;
		}
		
		FieldType calcElectronicEnergy(const std::vector<FieldType>& eigs) const
		{
			FieldType sum = 0;
			for (size_t i=0;i<eigs.size();i++) {
				sum += eigs[i]*utils::fermi((eigs[i]-engineParams_.mu)*engineParams_.beta);
			}
			return sum;
				
		}
		
		FieldType calcSuperExchange(const DynVarsType& dynVars) const
		{
			FieldType sum = 0;
			for (size_t i=0;i<geometry_.volume();i++) {
				for (size_t k = 0; k<geometry_.z(1); k++){
					size_t j=geometry_.neighbor(i,k).first;
					FieldType t1=dynVars.theta[i];
					FieldType t2=dynVars.theta[j];
					FieldType p1=dynVars.phi[i];
					FieldType p2=dynVars.phi[j];
					FieldType tmp = cos(t1)*cos(t2)+sin(t1)*sin(t2)*(cos(p1)*cos(p2)+sin(p1)*sin(p2));
					sum += mp_.jafNn*tmp;
				}
			}
			
			return sum*0.5;
		}
		
		FieldType calcMag(const DynVarsType& dynVars) const
		{
			std::vector<FieldType> mag(3);
			
			for (size_t i=0;i<geometry_.volume();i++) {
				mag[0] += sin(dynVars.theta[i])*cos(dynVars.phi[i]);
				mag[1] += sin(dynVars.theta[i])*sin(dynVars.phi[i]);
				mag[2] += cos(dynVars.theta[i]);
			}
			return (mag[0]*mag[0]+mag[1]*mag[1]+mag[2]*mag[2]);
		}
		
		FieldType calcKinetic(const DynVarsType& dynVars,
				      const std::vector<FieldType>& eigs) const
		{
			FieldType sum = 0;
			//const psimag::Matrix<ComplexType>& matrix = matrix_;
// 			for (size_t lambda=0;lambda<hilbertSize_;lambda++) {
// 				FieldType tmp2=0.0;
// 				for (size_t i=0;i<geometry_.volume();i++) {
// 					for (size_t k=0;k<geometry.z(1);k++) {
// 						size_t j=geometry.neighbor(i,k).first;
// 						for (size_t spin=0;spin<2;spin++) {
// 							ComplexType tmp = conj(matrix(i+spin*ether.linSize,lambda))*matrix(j+spin*ether.linSize,lambda);
// 							tmp2 += mp_.hoppings[orb1+spin*nbands_+dir*nbands_*nbands_]*real(tmp);
// 						}
// 					}
// 				}
// 				sum += tmp2 * utils::fermi(engineParams_.beta*(eigs[lambda]-engineParams_.mu));
// 			}
			return sum;
		}

		
		const EngineParamsType& engineParams_;
		const ParametersModelType& mp_;
		const GeometryType& geometry_;
		size_t hilbertSize_;
		MatrixType matrix_;
		ProgressIndicatorType progress_;
		RandomNumberGeneratorType rng_;
		
	}; // PnictidesTwoOrbitals
	
	template<typename FieldType,typename EngineParamsType,typename ParametersModelType,typename GeometryType>
	std::ostream& operator<<(std::ostream& os,const PnictidesTwoOrbitals<FieldType,EngineParamsType,ParametersModelType,GeometryType>& model)
	{
		os<<"ModelParameters\n";
		os<<model.mp_;
		return os;
	}
} // namespace Spf

/*@}*/
#endif // PNICTIDES_2ORB_H
