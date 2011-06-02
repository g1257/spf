
/** \ingroup SPF */
/*@{*/

/*! \file PnictidesMultiOrbitals
 *
 *  PnictidesMultiOrbitals model
 *
 */
#ifndef PNICTIDES_MULTI_ORBS_H
#define PNICTIDES_MULTI_ORBS_H
#include "PnictidesMultiOrbitalsFields.h"
#include "Random48.h"
#include "ProgressIndicator.h"
#include "Adjustments.h"
#include "SpinOperations.h"
#include "ModelBase.h"
#include "ThreeOrbitalTerms.h"
#include "ObservablesStored.h"
#include "Conductance.h"

namespace Spf {
	template<
		typename EngineParamsType,
		typename ParametersModelType_,
		typename GeometryType>
	class PnictidesMultiOrbitals : public ModelBase<Spin<typename EngineParamsType::FieldType>,EngineParamsType,ParametersModelType_,GeometryType> {
		
		typedef typename EngineParamsType::FieldType FieldType;
		typedef std::complex<FieldType> ComplexType;
		typedef PsimagLite::Matrix<ComplexType> MatrixType;
		//typedef RandomNumberGenerator<FieldType> RandomNumberGeneratorType;
		typedef typename GeometryType::PairType PairType;
		typedef Dmrg::ProgressIndicator ProgressIndicatorType;
		typedef Adjustments<EngineParamsType> AdjustmentsType;

		public:

		typedef ParametersModelType_ ParametersModelType;
		static const size_t norb_ = ParametersModelType::numberOfOrbitals;
		typedef PnictidesTwoOrbitalsFields<FieldType,GeometryType> DynVarsType;
		typedef typename DynVarsType::SpinType SpinType;
		typedef typename DynVarsType::SpinOperationsType SpinOperationsType;
		typedef ThreeOrbitalTerms<norb_,MatrixType,ParametersModelType,
				GeometryType> ThreeOrbitalTermsType;
		typedef ObservablesStored<SpinOperationsType,ComplexType,
				ParametersModelType> ObservablesStoredType;
		
		enum {OLDFIELDS,NEWFIELDS};
		enum {SPIN_UP,SPIN_DOWN};
		
		PnictidesMultiOrbitals(
				const EngineParamsType& engineParams,
				const ParametersModelType& mp,
				const GeometryType& geometry) :
			engineParams_(engineParams),mp_(mp),geometry_(geometry),dynVars_(geometry.volume(),engineParams.dynvarsfile),
				      hilbertSize_(2*mp_.numberOfOrbitals*geometry.volume()),
				      adjustments_(engineParams),progress_("PnictidesTwoOrbitals",0),
					spinOperations_(geometry,engineParams.mcWindow),
					threeOrbitalTerms_(mp,geometry),
					observablesStored_(spinOperations_,geometry,mp_,2*mp_.numberOfOrbitals)
		{
		}
		
		DynVarsType& dynVars() { return dynVars_; }
		
		size_t totalFlips() const { return geometry_.volume(); }
		
		SpinOperationsType& ops(SpinOperationsType*) { return spinOperations_; }
		
		size_t hilbertSize() const { return hilbertSize_; }
		
		FieldType deltaDirect(size_t i) const 
		{
			FieldType x = spinOperations_.deltaDirect(i,mp_.jafNn,mp_.jafNnn);
			x += spinOperations_.deltaMagneticField(i,mp_.magneticField);
			return x;
		}
		
		void set(typename DynVarsType::SpinType& dynVars) { spinOperations_.set(dynVars); }
		
		template<typename GreenFunctionType>
		void doMeasurements(GreenFunctionType& greenFunction,size_t iter,std::ostream& fout)
		{
			const SpinType& dynVars = dynVars_.getField((SpinType*)0);
			
			std::string s = "iter=" + ttos(iter); 
			progress_.printline(s,fout);
				
			FieldType temp=greenFunction.calcNumber();
			s ="Number_Of_Electrons="+ttos(temp);
			progress_.printline(s,fout);
			
			//s = "rankGlobal=";
			
			temp=greenFunction.calcElectronicEnergy();
			s="Electronic Energy="+ttos(temp);
			progress_.printline(s,fout);
			
			FieldType temp2=spinOperations_.calcSuperExchange(dynVars,mp_.jafNn);
			s="Superexchange="+ttos(temp2);
			progress_.printline(s,fout);
			
			temp += temp2;
			if (mp_.jafNnn!=0) {
				temp2=spinOperations_.directExchange2(dynVars,mp_.jafNnn);
				s="Superexchange2="+ttos(temp2);
				progress_.printline(s,fout);
				temp += temp2;
			}

			// total energy = electronic energy + superexchange + phonon energy
			s="TotalEnergy="+ttos(temp);
			progress_.printline(s,fout);
				
			//s="Action=";
			
			//s="Number_Of_Holes=";
			
			adjustments_.print(fout);
			
			std::vector<FieldType> magVector(3,0);
			spinOperations_.calcMagVector(magVector,dynVars);
			s="ClassicalMagnetizationSquared="+ttos(magVector*magVector);
			progress_.printline(s,fout);

			std::vector<ComplexType> electronSpinVector(3,0);
			greenFunction.electronSpin(electronSpinVector,mp_.numberOfOrbitals,dynVars.size);
			std::vector<ComplexType> combinedVector(3,0);
			combinedVector =  electronSpinVector + magVector;
			s="CombinedMagnetizationSquared="+ttos(std::real(combinedVector*combinedVector));
			progress_.printline(s,fout);
			
			for (size_t i = 0;i<combinedVector.size();i++) {
				s="CombinedMagnetization"+ttos(i)+"="+ttos(combinedVector[i]);
				progress_.printline(s,fout);
			}

			PsimagLite::Matrix<FieldType> v
				(greenFunction.hilbertSize(),greenFunction.hilbertSize());
			calcVelocitySquared(greenFunction,v);
			typedef Conductance<EngineParamsType,GreenFunctionType> ConductanceType;
			ConductanceType conductance(engineParams_,greenFunction);
			s = "Conductance=" + ttos(conductance(v));
			progress_.printline(s,fout);

// 			temp=calcKinetic(dynVars_,eigs);
// 			s ="KineticEnergy="+ttos(temp);
// 			progress_.printline(s,fout);
			
			observablesStored_(dynVars,greenFunction);
		} // doMeasurements
		
		void createHamiltonian(MatrixType& matrix,size_t oldOrNewDynVars)
		{
			const SpinType& dynVars = dynVars_.getField((SpinType*)0);
			if (oldOrNewDynVars==NEWFIELDS) createHamiltonian(spinOperations_.dynVars2(),matrix);
			else createHamiltonian(dynVars,matrix);
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
			spinOperations_.accept(i);
		}
		
		FieldType integrationMeasure(size_t i)
		{
			return spinOperations_.sineUpdate(i);
		}
		
		void finalize(std::ostream& fout)
		{
			observablesStored_.finalize(fout);	
		}
		
		template<typename EngineParamsType2,typename ParametersModelType2,typename GeometryType2>
		friend std::ostream& operator<<(std::ostream& os,
				const PnictidesMultiOrbitals<EngineParamsType2,ParametersModelType2,GeometryType2>& model);
		
		private:
		
		void createHamiltonian(const typename DynVarsType::SpinType& dynVars,MatrixType& matrix) const
		{
			size_t volume = geometry_.volume();
			size_t norb = mp_.numberOfOrbitals;
			size_t dof = norb * 2; // the 2 comes because of the spin
			std::vector<ComplexType> jmatrix(2*2);
			
			for (size_t gamma1=0;gamma1<matrix.n_row();gamma1++) 
				for (size_t p = 0; p < matrix.n_col(); p++) 
					matrix(gamma1,p)=0;
			
			for (size_t p = 0; p < volume; p++) {
				auxCreateJmatrix(jmatrix,dynVars,p);
				for (size_t gamma1=0;gamma1<dof;gamma1++) {		
					
		
					size_t spin1 = size_t(gamma1/norb);
					size_t orb1 = gamma1 % norb;
					//! Term B (n_iup - n_idown)
					FieldType magField = (spin1==SPIN_UP) ? mp_.magneticField : -mp_.magneticField;
					matrix(p+gamma1*volume,p+gamma1*volume) =
						real(jmatrix[spin1+2*spin1]) + mp_.potentialV[p] + magField;
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
						matrix(p + gamma2*volume,p+gamma1*volume) =
								conj(matrix(p+gamma1*volume,p + gamma2*volume));
					}
				}
			}
			threeOrbitalTerms_(matrix);
		}

		void auxCreateJmatrix(std::vector<ComplexType>& jmatrix,const
				typename DynVarsType::SpinType& dynVars,size_t site) const
		{
			if (!mp_.modulus[site]) {
				for (size_t i=0;i<jmatrix.size();i++) jmatrix[i] = 0;
				return;
			}
			
			jmatrix[0]=cos(dynVars.theta[site]);
		
			jmatrix[1]=ComplexType(sin(dynVars.theta[site])*cos(dynVars.phi[site]),
				sin(dynVars.theta[site])*sin(dynVars.phi[site]));
		
			jmatrix[2]=conj(jmatrix[1]);
		
			jmatrix[3]= -cos(dynVars.theta[site]);
		
			for (size_t i=0;i<jmatrix.size();i++) jmatrix[i] *= mp_.J;
		}
		
		FieldType calcKinetic(const DynVarsType& dynVars,
				      const std::vector<FieldType>& eigs) const
		{
			FieldType sum = 0;
			//const PsimagLite::Matrix<ComplexType>& matrix = matrix_;
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
// 				sum += tmp2 * fermi(engineParams_.beta*(eigs[lambda]-engineParams_.mu));
// 			}
			return sum;
		}

		//! Assuming a constant hopping for all spin and orbitals
		//! FIXME: Replace with actual velocity for this model
		template<typename GreenFunctionType>
		void calcVelocitySquared(const GreenFunctionType& gf,
				PsimagLite::Matrix<FieldType>& v) const
		{
			size_t ly = geometry_.length();
			size_t norb = mp_.numberOfOrbitals;
			size_t volume = geometry_.volume();
			size_t dof = 2 * norb; // 2 is for the spin

			for (size_t a=0;a<v.n_row();a++) {
				for (size_t b=0;b<v.n_col();b++) {
					ComplexType sum = 0;
					for (size_t y=0;y<ly;y++) {
						for (size_t gamma=0;gamma<dof;gamma++) {
							size_t i = geometry_.coorToIndex(0,y); // x=0;
							size_t j = geometry_.coorToIndex(1,y); // x=1;
							i += gamma*volume; // add spin and orb.
							j += gamma*volume; // add spin and orb.
							sum += velocity(gf,i,j,a,b);
						}
					}
					v(a,b) = std::real(std::conj(sum)*sum);
				}
			}
		}

		//! Assuming a constant hopping for all spin and orbitals
		//! FIXME: Replace with actual velocity for this model
		template<typename GreenFunctionType>
		ComplexType velocity(const GreenFunctionType& gf,
				size_t i,
				size_t j,
				size_t a,
				size_t b) const
		{
			return std::conj(gf.matrix(i,a)) * gf.matrix(j,b) -
					std::conj(gf.matrix(j,a)) * gf.matrix(i,a);
		}

		const EngineParamsType& engineParams_;
		const ParametersModelType& mp_;
		const GeometryType& geometry_;
		DynVarsType dynVars_;
		size_t hilbertSize_;
		AdjustmentsType adjustments_;
		ProgressIndicatorType progress_;
		SpinOperationsType spinOperations_;
		ThreeOrbitalTermsType threeOrbitalTerms_;
		ObservablesStoredType observablesStored_;

	}; // PnictidesMultiOrbitals

	template<typename EngineParamsType,typename ParametersModelType,typename GeometryType>
	std::ostream& operator<<(std::ostream& os,const PnictidesMultiOrbitals<EngineParamsType,ParametersModelType,GeometryType>& model)
	{
		os<<"ModelParameters\n";
		os<<model.mp_;
		return os;
	}
} // namespace Spf

/*@}*/
#endif // PNICTIDES_MULTI_ORBS_H
