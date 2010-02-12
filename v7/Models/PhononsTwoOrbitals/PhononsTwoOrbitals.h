
/** \ingroup SPF */
/*@{*/

/*! \file PhononsTwoOrbitals.h
 *
 *  PhononsTwoOrbitals model (JH infinity)
 *
 */
#ifndef PHONONS_2ORB_H
#define PHONONS_2ORB_H
#include "Utils.h"
#include "Spin.h"
#include "RandomNumberGenerator.h"
#include "ProgressIndicator.h"
#include "Adjustments.h"
#include "SpinOperations.h"
#include "PhononOperations.h"
#include "MonteCarlo.h"
#include "ModelBase.h"

namespace Spf {
	template<typename EngineParamsType,typename ParametersModelType_,typename GeometryType>
	class PhononsTwoOrbitals : public ModelBase
			<Spin<typename EngineParamsType::FieldType>,EngineParamsType,ParametersModelType_,GeometryType> {
		
		typedef typename EngineParamsType::FieldType FieldType;
		typedef std::complex<FieldType> ComplexType;
		typedef psimag::Matrix<ComplexType> MatrixType;
		typedef typename GeometryType::PairType PairType;
		typedef Dmrg::ProgressIndicator ProgressIndicatorType;
		typedef Adjustments<EngineParamsType> AdjustmentsType;
		typedef PhononsTwoOrbitals<EngineParamsType,ParametersModelType_,GeometryType> ThisType;
		
		static const size_t nbands_ = 2;
		
		public:
		typedef ParametersModelType_ ParametersModelType;
		typedef Spin<FieldType> DynVarsType;
		typedef ClassicalSpinOperations<GeometryType,DynVarsType> ClassicalSpinOperationsType;
		typedef PhononOperations<GeometryType,DynVarsType> PhononOperationsType;
		
		enum {OLDFIELDS,NEWFIELDS};
		
		PhononsTwoOrbitals(const EngineParamsType& engineParams,const ParametersModelType& mp,const GeometryType& geometry) :
			engineParams_(engineParams),mp_(mp),geometry_(geometry),dynVars_(geometry.volume(),engineParams.dynvarsfile),
				      hilbertSize_(nbands_*geometry_.volume()), // there's no spin here
				      adjustments_(engineParams),progress_("PhononsTwoOrbitals",0),
					classicalSpinOperations_(geometry_,engineParams_.mcWindow)
		{
		}
		
		DynVarsType& dynVars() { return dynVars_; }
		
		size_t totalFlips() const { return geometry_.volume(); }
		
		//size_t dof() const { return 2; }
		
		size_t hilbertSize() const { return hilbertSize_; }
		
		FieldType deltaDirect(size_t i) const 
		{
			return classicalSpinOperations_.deltaDirect(i,mp_.jaf,0);
		}
		
		void set(DynVarsType& dynVars) { classicalSpinOperations_.set(dynVars_); }
		
		template<typename RandomNumberGeneratorType>
		void propose(size_t i,RandomNumberGeneratorType& rng) { classicalSpinOperations_.propose(i,rng); }
				
		template<typename GreenFunctionType>
		void doMeasurements(const DynVarsType& dynVars,GreenFunctionType& greenFunction,size_t iter,std::ostream& fout)
		{
			std::string s = "iter=" + utils::ttos(iter); 
			progress_.printline(s,fout);
				
			FieldType temp=calcNumber(greenFunction);
			s ="Number_Of_Electrons="+utils::ttos(temp);
			progress_.printline(s,fout);
			
			//s = "rankGlobal=";
			
			temp=calcElectronicEnergy(greenFunction);
			s="Electronic Energy="+utils::ttos(temp);
			progress_.printline(s,fout);
			
			FieldType temp2=calcSuperExchange(dynVars_);
			s="Superexchange="+utils::ttos(temp2);
			progress_.printline(s,fout);
			
			temp += temp2;

			// total energy = electronic energy + superexchange + phonon energy
			s="TotalEnergy="+utils::ttos(temp);
			progress_.printline(s,fout);
				
			//s="Action=";
			
			//s="Number_Of_Holes=";
			
			adjustments_.print(fout);
			
			temp = calcMag(dynVars_);
			s="Mag2="+utils::ttos(temp);
			progress_.printline(s,fout);
			
// 			temp=calcKinetic(dynVars_,eigs);
// 			s ="KineticEnergy="+utils::ttos(temp);
// 			progress_.printline(s,fout);
			
			//storedObservables_.doThem();
		} // doMeasurements
		
		void createHamiltonian(psimag::Matrix<ComplexType>& matrix,size_t oldOrNewDynVars) const
		{
			 if (oldOrNewDynVars==NEWFIELDS) createHamiltonian(classicalSpinOperations_.dynVars2(),matrix);
			 else createHamiltonian(dynVars_,matrix);
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
		
		template<typename EngineParamsType2,typename ParametersModelType2,typename GeometryType2>
		friend std::ostream& operator<<(std::ostream& os,
				const PhononsTwoOrbitals<EngineParamsType2,ParametersModelType2,GeometryType2>& model);
		
		private:
		
		void createHamiltonian(const DynVarsType& dynVars,MatrixType& matrix) const
		{
			size_t volume = geometry_.volume();
			size_t norb = nbands_;
			size_t dof = norb; // no spin here
			
			for (size_t gamma1=0;gamma1<matrix.n_row();gamma1++) 
				for (size_t p = 0; p < matrix.n_col(); p++) 
					matrix(gamma1,p)=0;

			for (size_t p = 0; p < volume; p++) {
				FieldType phonon_q1=phononOperations_.calcPhonon(p,dynVars,0);
				FieldType phonon_q2=phononOperations_.calcPhonon(p,dynVars,1);
				FieldTyp phonon_q3=phononOperations_.calcPhonon(p,dynVars,2);	
				matrix(p,p) = mp_.phononSpinCouplings[0]*phonon_q1+
						mp_.phononSpinCouplings[2]*phonon_q3+
						mp_.potential[p];
				matrix(p+volume,p+volume) = -mp_.phononSpinCouplings[2]*phonon_q3+
						mp_.phononSpinCouplings[0]*phonon_q1+mp_.potential[p];
				matrix(p,p+volume) = (mp_.phononSpinCouplings[1]*phonon_q2);
				matrix(p+volume,p) = conj(matrix(p,p+volume));
				
				for (size_t j = 0; j < geometry.z(1); j++) {	/* hopping elements, n-n only */
					PairType tmpPair = geometry_.neighbor(p,j);
					size_t col = tmpPair.first;
					size_t dir = tmpPair.second; // int(j/2);
					
					FieldType tmp=cos(0.5*dynVars.theta[p])*cos(0.5*dynVars.theta[col]);
					FieldType tmp2=sin(0.5*dynVars.theta[p])*sin(0.5*dynVars.theta[col]);
					FieldType S_ij=tpem_t(tmp+tmp2*cos(dynVars.phi[p]-dynVars.phi[col]),
						-tmp2*sin(dynVars.phi[p]-dynVars.phi[col]));
					
					matrix(p,col) = -mp_.bandHoppings[0+0*2+dir*4] * S_ij;
					matrix(col,p) = conj(hopping * S_ij);
					
					matrix(p, col+volume)= -mp_.bandHoppings[0+1*2+dir*4] * S_ij;
					matrix(col+volume,p)=conj(matrix(p,col+volume));
					
					matrix(p+volume,col) =  -mp_.bandHoppings[1+0*2+dir*4] * S_ij;
					matrix(col,p+volume) =  conj(matrix(p+volume,col));
					
					matrix(p+volume,col+volume) =  -ether.bandHoppings[1+1*2+dir*4] * S_ij;
					matrix(col+volume,p+volume) = conj(matrix(p+volume,col+volume));
				}
			}
		}

		template<typename GreenFunctionType>
		FieldType calcNumber(GreenFunctionType& greenFunction) const
		{
			FieldType sum=0;
			for (size_t i=0;i<hilbertSize_;i++) {
				sum += real(greenFunction(i,i));
			}
			return sum;
		}

		template<typename GreenFunctionType>
		FieldType calcElectronicEnergy(GreenFunctionType& greenFunction) const
		{
			FieldType sum = 0;
			for (size_t i=0;i<hilbertSize_;i++) {
				FieldType g = real(greenFunction(i,i));
				sum += g*(engineParams_.mu+log(1./g-1.)/engineParams_.beta);
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
		DynVarsType dynVars_;
		size_t hilbertSize_;
		AdjustmentsType adjustments_;
		ProgressIndicatorType progress_;
		//RandomNumberGeneratorType& rng_;
		ClassicalSpinOperationsType classicalSpinOperations_;
		PhononsType phononOperations_;
	}; // PhononsTwoOrbitals

	template<typename EngineParamsType,typename ParametersModelType,typename GeometryType>
	std::ostream& operator<<(std::ostream& os,const PhononsTwoOrbitals<EngineParamsType,ParametersModelType,GeometryType>& model)
	{
		os<<"ModelParameters\n";
		os<<model.mp_;
		return os;
	}
} // namespace Spf

/*@}*/
#endif // PHONONS_2ORB_H
