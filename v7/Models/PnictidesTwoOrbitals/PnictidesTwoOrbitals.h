
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
#include "Spin.h"
#include "RandomNumberGenerator.h"
#include "ProgressIndicator.h"
#include "Adjustments.h"
#include "SpinOperations.h"
#include "MonteCarlo.h"

namespace Spf {
	template<typename FieldType,typename EngineParamsType,typename ParametersModelType,typename GeometryType>
	class PnictidesTwoOrbitals {
		
		typedef std::complex<FieldType> ComplexType;
		typedef psimag::Matrix<ComplexType> MatrixType;
		typedef RandomNumberGenerator<FieldType> RandomNumberGeneratorType;
		typedef typename GeometryType::PairType PairType;
		typedef Dmrg::ProgressIndicator ProgressIndicatorType;
		typedef Adjustments<EngineParamsType> AdjustmentsType;
		typedef PnictidesTwoOrbitals<FieldType,EngineParamsType,ParametersModelType,GeometryType> ThisType;
		
		static const size_t nbands_ = 2;
		
		public:
			
		typedef DynVars<FieldType> DynVarsType;
		typedef ClassicalSpinOperations<GeometryType,RandomNumberGeneratorType,DynVarsType> ClassicalSpinOperationsType;
		typedef MonteCarlo<EngineParamsType,ThisType,DynVarsType,RandomNumberGeneratorType> MonteCarloType;
		
		PnictidesTwoOrbitals(const EngineParamsType& engineParams,const ParametersModelType& mp,const GeometryType& geometry) :
			engineParams_(engineParams),mp_(mp),geometry_(geometry),hilbertSize_(2*nbands_*geometry_.volume()),
				      matrix_(hilbertSize_,hilbertSize_),adjustments_(engineParams),progress_("PnictidesTwoOrbitals",0),
					rng_(),classicalSpinOperations_(geometry_,rng_,engineParams_.mcWindow),
					monteCarlo_(engineParams,*this,rng_)
		{
		}
		
		size_t totalFlips() const { return geometry_.volume(); }
		
		FieldType deltaDirect(size_t i) const 
		{
			return classicalSpinOperations_.deltaDirect(i,mp_.jafNn,mp_.jafNnn);
		}
		
		void set(DynVarsType& dynVars) { classicalSpinOperations_.set(dynVars); }
		
		void propose(size_t i) { classicalSpinOperations_.propose(i); }
		

		size_t doMonteCarlo(DynVarsType& dynVars, size_t iter)
		{
			return monteCarlo_(dynVars,iter);
			
		}
				
		void doMeasurements(DynVarsType& dynVars, size_t iter,std::ostream& fout)
		{
			std::string s = "iter=" + utils::ttos(iter); 
			
			progress_.printline(s,fout);
				
			std::vector<FieldType> eigs;
			fillAndDiag(eigs,dynVars,'V');
				
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
				temp2=classicalSpinOperations_.directExchange2(dynVars,mp_.jafNnn);
				s="Superexchange2="+utils::ttos(temp2);
				progress_.printline(s,fout);
				temp += temp2;
			}

			// total energy = electronic energy + superexchange + phonon energy
			s="TotalEnergy="+utils::ttos(temp);
			progress_.printline(s,fout);
				
			//s="Action=";
			
			//s="Number_Of_Holes=";
			
			adjustments_.print(fout);
			
			temp = calcMag(dynVars);
			s="Mag2="+utils::ttos(temp);
			progress_.printline(s,fout);
			
			temp=calcKinetic(dynVars,eigs);
			s ="KineticEnergy="+utils::ttos(temp);
			progress_.printline(s,fout);
			
			//storedObservables_.doThem();
		} // doMeasurements
		
		void fillAndDiag(std::vector<FieldType> &eig)
		{
			fillAndDiag(eig,classicalSpinOperations_.dynVars2());	
		}
		
		void fillAndDiag(std::vector<FieldType> &eig,const DynVarsType& dynVars,char jobz='N')
		{
			if (jobz=='v') jobz='V';
			
			createHamiltonian(dynVars,matrix_);
			utils::diag(matrix_,eig,jobz);

			if (jobz!='V') sort(eig.begin(), eig.end(), std::less<FieldType>());
		}
		
		void adjustChemPot(const std::vector<FieldType>& eigs)
		{
			if (engineParams_.carriers==0) return;
			try {
				engineParams_.mu = adjustments_.adjChemPot(eigs);
			} catch (std::exception& e) {
				std::cerr<<e.what()<<"\n";
			}
				
		}
		
		void accept(size_t i) 
		{
			return classicalSpinOperations_.accept(i);
		}
		
		FieldType integrationMeasure(size_t i)
		{
			return classicalSpinOperations_.sineUpdate(i);
		}
		
		template<typename FieldType2,typename EngineParamsType2,typename ParametersModelType2,typename GeometryType2>
		friend std::ostream& operator<<(std::ostream& os,
				const PnictidesTwoOrbitals<FieldType2,EngineParamsType2,ParametersModelType2,GeometryType2>& model);
		
		private:
		
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
		AdjustmentsType adjustments_;
		ProgressIndicatorType progress_;
		RandomNumberGeneratorType rng_;
		ClassicalSpinOperationsType classicalSpinOperations_;
		MonteCarloType monteCarlo_;
		
		
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
