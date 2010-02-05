
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

namespace Spf {
	class PnictidesTwoOrbitals {
		public:
		size_t doMonteCarlo(size_t iter)
		{
			size_t acc = 0;
			for (size_t i=0;i<n;i++) {
				
				dynVars2.theta=dynVars.theta;
				dynVars2.phi=dynVars.phi;
					
				r_newSpins(dynVars.theta[i],dynVars.phi[i],thetaNew,phiNew);
				
				dynVars2.theta[i]=thetaNew;
				dynVars2.phi[i]=phiNew;
				
				dsDirect = dSDirect(dynVars,dynVars2,i);
				dsDirect += directExchange2(dynVars2)-directExchange2(dynVars);
				
				oldmu=aux.varMu;	
				diag(eigNewAllBands,geometry,dynVars2,ether,aux);
				if (ether.carriers>0) adjChemPot(eigNewAllBands,ether,aux);
				sineupdate= sin(dynVars.theta[i]);
				if (dynVars.theta[i]!=0) {
					sineupdate = sin(dynVars2.theta[i])/sineupdate;
				} else {
					sineupdate = 1.0;
				}
				flag=doMetropolis(eigNewAllBands,aux.eigAllBands,dsDirect,sineupdate);
					
				if (flag && (mp_.mcflag & 1)) { // Accepted
					dynVars.theta[i]=thetaNew;
					dynVars.phi[i]=phiNew;
					acc++;
				} else { // not accepted
						//aux.varMu=oldmu;
				}
			} // lattice sweep
			return acc;
		}
		private:
			
		void r_newSpins(FieldType thetaOld, FieldType phiOld, FieldType &thetaNew,FieldType &phiNew)
		{
			if (isingSpins_) {
				if (thetaOld==0) thetaNew=M_PI; 
				else thetaNew=0;
				phiNew=0;
			} else {
				if (mp_.mcwindow<0) {
					thetaNew = 2*rng_()-1;
					phiNew = 2*M_PI*rng_();
					thetaNew = acos(thetaNew);
				} else {
					thetaNew=2*rng_()- 1;
					if (thetaNew < -1) thetaNew= 0;
					if (thetaNew > 1) thetaNew = 0;		
					phiNew=phiOld+2*M_PI*(rng_()- 0.5)*mp_.mcwindow;
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
			size_t n = geometry.volume();
			FieldType dS = 0;
			
			for (size_t i=0;i<n;i++) {
				for (size_t k = 0; k<geometry_.z(i,2); k++){
					size_t j=geometry_.neighbor(i,k,2); /**next nearest neighbor */
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
		
		FieldType dSDirect(const DynVars& dynVars,const DynVars& dynVars2, size_t i)
		{
			FieldType dS = 0;
				
			for (size_t k = 0; k<geometry.z(i); k++){
				size_t j=geometry.neighbor(i,k);
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
		
		void diag(std::vector<FieldType> &eig,const DynVars& dynVars,char jobz)
		{
			if (jobz=='v') jobz='V';
	
			createHamiltonian(dynVars,matrix);
			utils::diag(matrix_,eig,jobz);

			if (jobz!='V') sort(eig.begin(), eig.end(), less<double>());
		}
		
		bool doMetropolis(const std::vector<FieldType>& eNew,const std::vector<FieldType>& eOld,
			FieldType dsDirect,FieldType sineupdate)
		{
			FieldType mu=params_.mu;
			FieldType beta = params_.beta;
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
		
		void createHamiltonian (const DynVars& dynVars,MatrixType& matrix)
		{
			size_t volume = geometry_.volume();
			size_t norb = mp_.numberOfOrbitals;
			size_t dof = norb * 2; // the 2 comes because of the spin
			std::vector<FieldType> jmatrix(dof*dof);
			
			for (size_t gamma1=0;gamma1<matrix.n_row();gamma1++) 
				for (size_t p = 0; p < matrix.n_col(); p++) 
					matrix(gamma1,p)=0;
			
			for (size_t gamma1=0;gamma1<dof;gamma1++) {
				for (size_t p = 0; p < volume; p++) {
					
					auxCreateJmatrix(jmatrix,dynVars,ether,p);
		
					size_t spin1 = size_t(gamma1/2);
					size_t orb1 = gamma1 % 2;
					matrix(p+gamma1*volume,p+gamma1*volume) = real(jmatrix[spin1+2*spin1]) + ether.potential[p];
						
					for (size_t j = 0; j <  geometry.z(p); j++) {	
						if (j%2!=0) continue;	
						PairType tmpPair = geometry.neighbor(p,j);
						size_t k = tmpPair.first;
						size_t dir = tmpPair.second; // int(j/2);
						for (size_t orb2=0;orb2<norb;orb2++) {
							size_t gamma2 = orb2+spin1*norb; // spin2 == spin1 here
							matrix(p+gamma1*volume,k+gamma2*volume)=
								ether.bandHoppings[orb1+orb2*norb+norb*norb*dir];
							matrix(k+gamma2*volume,p+gamma1*volume) = conj(matrix(p+gamma1*volume,k+gamma2*volume));
						}
					}
					//if (geometry.z(p,2)!=4 || geometry.z(p)!=4) throw std::runtime_error("neighbours wrong\n");
					for (size_t j = 0; j <  geometry.z(p,2); j++) {
						if (j%2!=0) continue;	
						PairType tmpPair = geometry.neighbor(p,j,2);
						size_t k = tmpPair.first;
						size_t dir = tmpPair.second + 2; // int dir = int(j/2)+2;
						//std::cerr<<"Neigbors "<<p<<" "<<k<<"\n";
						for (size_t orb2=0;orb2<norb;orb2++) {
							size_t gamma2 = orb2+spin1*norb; // spin2 == spin1 here
							matrix(p+gamma1*volume,k+gamma2*volume)=
								ether.bandHoppings[orb1+orb2*norb+norb*norb*dir];
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
		
		void auxCreateJmatrix(std::vector<FieldType>& jmatrix,const DynVarsType& dynVars,size_t site)
		{
			
			jmatrix[0]=0.5*cos(dynVars.theta[site]);
		
			jmatrix[1]=MatType(sin(dynVars.theta[site])*cos(dynVars.phi[site]),
				sin(dynVars.theta[site])*sin(dynVars.phi[site]));     
		
			jmatrix[2]=conj(jmatrix[1]);
		
			jmatrix[3]= -cos(dynVars.theta[site]);
		
			for (size_t i=0;i<jmatrix.size();i++) jmatrix[i] *= ether.J[0]; 
		}
		
	}; // PnictidesTwoOrbitals
} // namespace Spf

/*@}*/
#endif // PNICTIDES_2ORB_H
