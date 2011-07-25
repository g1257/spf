
/** \ingroup SPF */
/*@{*/

/*! \file DmsMultiOrbital.h
 *
 *  DMS multiorbital model
 *
 */
#ifndef DMS_MULIORBITAL_H
#define DMS_MULIORBITAL_H
#include "DmsMultiOrbitalFields.h"
#include "Random48.h"
#include "ProgressIndicator.h"
#include "Adjustments.h"
#include "SpinOperations.h"
#include "ModelBase.h"
#include "ObservablesStored.h"

namespace Spf {
	template<
		typename EngineParamsType,
		typename ParametersModelType_,
		typename GeometryType,
		typename ConcurrencyType>
	class DmsMultiOrbital : public ModelBase<Spin<
	     typename EngineParamsType::FieldType>,
		 EngineParamsType,ParametersModelType_,GeometryType,ConcurrencyType> {
		
		typedef typename EngineParamsType::FieldType RealType;
		typedef std::complex<RealType> ComplexType;
		typedef PsimagLite::Matrix<ComplexType> MatrixType;
		//typedef RandomNumberGenerator<RealType> RandomNumberGeneratorType;
		typedef typename GeometryType::PairType PairType;
		typedef Dmrg::ProgressIndicator ProgressIndicatorType;
		typedef Adjustments<EngineParamsType> AdjustmentsType;
		typedef DmsMultiOrbital<EngineParamsType,ParametersModelType_,
		                        GeometryType,ConcurrencyType> ThisType;

	public:
		typedef ParametersModelType_ ParametersModelType;
		typedef DmsMultiOrbitalFields<RealType,GeometryType> DynVarsType;
		typedef typename DynVarsType::SpinType SpinType;
		typedef typename DynVarsType::SpinOperationsType SpinOperationsType;
		typedef ObservablesStored<SpinOperationsType,ComplexType,
		                          ParametersModelType,
		                          EngineParamsType,
		                          ConcurrencyType> ObservablesStoredType;
		static const size_t ORBITALS = ObservablesStoredType::ORBITALS;
		
		enum {OLDFIELDS,NEWFIELDS};
		
		DmsMultiOrbital(const EngineParamsType& engineParams,
		                const ParametersModelType& mp,
		                const GeometryType& geometry,
		                ConcurrencyType& concurrency)
		: engineParams_(engineParams),mp_(mp),
		 geometry_(geometry),
		 dynVars_(geometry.volume(),engineParams.dynvarsfile),
		 hilbertSize_(2*ORBITALS*geometry.volume()),
		 adjustments_(engineParams),progress_("PnictidesTwoOrbitals",0),
		 spinOperations_(geometry,engineParams.mcWindow),
		 observablesStored_(spinOperations_,geometry,mp_,engineParams_,concurrency)
		{
		}
		
		DynVarsType& dynVars() { return dynVars_; }
		
		size_t totalFlips() const { return geometry_.volume(); }
		
		SpinOperationsType& ops(SpinOperationsType*) { return spinOperations_; }
		
		size_t hilbertSize() const { return hilbertSize_; }
		
		RealType deltaDirect(size_t i) const
		{
			// don't use line below, unless geometry has nnn:
			//return spinOperations_.deltaDirect(i,mp_.jafNn,mp_.jafNnn);
			return 0.0;
		}
		
		void set(typename DynVarsType::SpinType& dynVars)
		{
			spinOperations_.set(dynVars);
		}
		
		template<typename GreenFunctionType>
		void doMeasurements(
				GreenFunctionType& greenFunction,
				size_t iter,
				std::ostream& fout)
		{
			const SpinType& dynVars = dynVars_.getField((SpinType*)0);
			
			std::string s = "iter=" + ttos(iter); 
			progress_.printline(s,fout);
				
			RealType temp=greenFunction.calcNumber();
			s ="Number_Of_Electrons="+ttos(temp);
			progress_.printline(s,fout);
			
			//s = "rankGlobal=";
			
			temp=greenFunction.calcElectronicEnergy();
			s="Electronic Energy="+ttos(temp);
			progress_.printline(s,fout);
			
			RealType temp2=spinOperations_.calcSuperExchange(dynVars,mp_.jafNn);
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
			
			temp = spinOperations_.calcMag(dynVars,mp_.modulus);
			s="Mag2="+ttos(temp);
			progress_.printline(s,fout);

			
// 			temp=calcKinetic(dynVars_,eigs);
// 			s ="KineticEnergy="+ttos(temp);
// 			progress_.printline(s,fout);
			
			observablesStored_(dynVars,greenFunction);
		} // doMeasurements
		
		void createHamiltonian(
				PsimagLite::Matrix<ComplexType>& matrix,
				size_t oldOrNewDynVars)
		{
			const SpinType& dynVars = dynVars_.getField((SpinType*)0);
			if (oldOrNewDynVars==NEWFIELDS)
				createHamiltonian(spinOperations_.dynVars2(),matrix);
			else
				createHamiltonian(dynVars,matrix);
		}
		
		void adjustChemPot(const std::vector<RealType>& eigs)
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
		
		RealType integrationMeasure(size_t i)
		{
			return spinOperations_.sineUpdate(i);
		}
		
		void finalize(std::ostream& fout)
		{
			observablesStored_.finalize(fout);	
		}
		
		template<
		         typename EngineParamsType2,
		         typename ParametersModelType2,
		         typename GeometryType2,
		         typename ConcurrencyType2>
		friend std::ostream& operator<<(
		           std::ostream& os,
		           const DmsMultiOrbital<EngineParamsType2,
		           ParametersModelType2,GeometryType2,ConcurrencyType2>& model);
		
	private:
		
		void createHamiltonian(
				const typename DynVarsType::SpinType& dynVars,
				MatrixType& matrix) const
		{
			size_t volume = geometry_.volume();
			size_t norb = ORBITALS;
			size_t dof = norb * 2; // the 2 comes because of the spin
			std::vector<ComplexType> jmatrix(dof*dof);
			std::vector<RealType> ymatrix(dof*dof);

			for (size_t gamma1=0;gamma1<matrix.n_row();gamma1++) 
				for (size_t p = 0; p < matrix.n_col(); p++) 
					matrix(gamma1,p)=0;

			for (size_t p = 0; p < volume; p++) {
				RealType modulus = mp_.modulus[p];

				auxCreateJmatrix(jmatrix,dynVars,p);
				auxCreateYmatrix(ymatrix,dynVars,p);
				for (size_t gamma1=0;gamma1<dof;gamma1++) {

					matrix(p+gamma1*volume,p+gamma1*volume) =
							real(jmatrix[gamma1+dof*gamma1])*modulus +
							mp_.potentialV[p] +
							ymatrix[gamma1+dof*gamma1];

					for (size_t j = 0; j <  geometry_.z(1); j++) {	
						//if (j%2!=0) continue;
						PairType tmpPair = geometry_.neighbor(p,j);
						size_t k = tmpPair.first;
						size_t dir = tmpPair.second; // int(j/2);
						for (size_t gamma2=0;gamma2<dof;gamma2++) {
							matrix(p+gamma1*volume,k+gamma2*volume)=
								mp_.hoppings[gamma1+gamma2*dof+dir*dof*dof];
							matrix(k+gamma2*volume,p+gamma1*volume) =
								conj(matrix(p+gamma1*volume,k+gamma2*volume));
						}
					}

					for (size_t gamma2=0;gamma2<dof;gamma2++) {
						matrix(p+gamma1*volume,p + gamma2*volume) =
							jmatrix[gamma1+dof*gamma2] * modulus;
						matrix(p + gamma2*volume,p+gamma1*volume) =
							conj(matrix(p+gamma1*volume,p + gamma2*volume));
					}
				}
			}
		}

		void auxCreateJmatrix(std::vector<ComplexType>& jmatrix,const
				typename DynVarsType::SpinType& dynVars,size_t site) const
		{
			jmatrix[0]=0.5*cos(dynVars.theta[site]);

			jmatrix[1]=0.0;

			jmatrix[2]=ComplexType(0.5*sin(dynVars.theta[site])*
					cos(dynVars.phi[site])/sqrt(3.0),
					-0.5*sin(dynVars.theta[site])*
					sin(dynVars.phi[site])/sqrt(3.0));

			jmatrix[3]=0.0;

			jmatrix[4]=ComplexType(sin(dynVars.theta[site])*
					cos(dynVars.phi[site])/sqrt(6.0),
					-sin(dynVars.theta[site])*
					sin(dynVars.phi[site])/sqrt(6.0));

			jmatrix[5]=0.0;
			jmatrix[6]=0.0;

			jmatrix[7]=-cos(dynVars.theta[site])/6.0;

			jmatrix[8]=ComplexType(sin(dynVars.theta[site])*
					cos(dynVars.phi[site])/3.0,
					sin(dynVars.theta[site])*
					sin(dynVars.phi[site])/3.0);

			jmatrix[9]=ComplexType(0.5*sin(dynVars.theta[site])*
					cos(dynVars.phi[site])/sqrt(3.0),
					-0.5*sin(dynVars.theta[site])*
					sin(dynVars.phi[site])/sqrt(3.0));

			jmatrix[10]=ComplexType(-sin(dynVars.theta[site])*
					cos(dynVars.phi[site])/(3.0*sqrt(2.0)),
					-sin(dynVars.theta[site])*
					sin(dynVars.phi[site])/(3.0*sqrt(2.0)));

			jmatrix[11]=-sqrt(2.0)*cos(dynVars.theta[site])/3.0;

			//jmatrix[11]=-2.0*cos(dynVars.theta[site])/sqrt(3.0);
			//         jmatrix[11]=0;
			jmatrix[12]= ComplexType(0.5*sin(dynVars.theta[site])*
					cos(dynVars.phi[site])/sqrt(3.0),
					0.5*sin(dynVars.theta[site])*
					sin(dynVars.phi[site])/sqrt(3.0));

			jmatrix[13]=ComplexType(sin(dynVars.theta[site])*
					cos(dynVars.phi[site])/3.0,
					-sin(dynVars.theta[site])*
					sin(dynVars.phi[site])/3.0);

			jmatrix[14]= -jmatrix[7];
			jmatrix[15]=0.0;
			jmatrix[16]= jmatrix[11];

			jmatrix[17]=ComplexType(sin(dynVars.theta[site])*
					cos(dynVars.phi[site])/(3.0*sqrt(2.0)),
					-sin(dynVars.theta[site])*
					sin(dynVars.phi[site])/(3.0*sqrt(2.0)));

			jmatrix[18]=0.0;

			jmatrix[19]=ComplexType(0.5*sin(dynVars.theta[site])*
					cos(dynVars.phi[site])/sqrt(3.0),
					0.5*sin(dynVars.theta[site])*
					sin(dynVars.phi[site])/sqrt(3.0));

			jmatrix[20]= 0.0;
			jmatrix[21]= -jmatrix[0];
			jmatrix[22]= 0.0;

			jmatrix[23]=ComplexType(-sin(dynVars.theta[site])*
					cos(dynVars.phi[site])/sqrt(6.0),
					-sin(dynVars.theta[site])*
					sin(dynVars.phi[site])/sqrt(6.0));

			jmatrix[24]=  ComplexType(sin(dynVars.theta[site])*
					cos(dynVars.phi[site])/sqrt(6.0),
					sin(dynVars.theta[site])*
					sin(dynVars.phi[site])/sqrt(6.0));

			jmatrix[25]=ComplexType(-sin(dynVars.theta[site])*
					cos(dynVars.phi[site])/(3.0*sqrt(2.0)),
					sin(dynVars.theta[site])*
					sin(dynVars.phi[site])/(3.0*sqrt(2.0)));

			jmatrix[26]=jmatrix[11];
			jmatrix[27]= 0.0;
			jmatrix[28]= jmatrix[7];

			jmatrix[29]= ComplexType(-sin(dynVars.theta[site])*
					cos(dynVars.phi[site])/6.0,
					sin(dynVars.theta[site])*
					sin(dynVars.phi[site])/6.0);

			jmatrix[30]= 0.0;
			jmatrix[31]= jmatrix[11];

			jmatrix[32]=ComplexType(sin(dynVars.theta[site])*
					cos(dynVars.phi[site])/(3.0*sqrt(2.0)),
					sin(dynVars.theta[site])*
					sin(dynVars.phi[site])/(3.0*sqrt(2.0)));

			jmatrix[33]=ComplexType(-sin(dynVars.theta[site])*
					cos(dynVars.phi[site])/sqrt(6.0),
					sin(dynVars.theta[site])*
					sin(dynVars.phi[site])/sqrt(6.0));

			jmatrix[34]= ComplexType(-sin(dynVars.theta[site])*
					cos(dynVars.phi[site])/6.0,
					-sin(dynVars.theta[site])*
					sin(dynVars.phi[site])/6.0);

			jmatrix[35]= jmatrix[14];

			size_t dof = ORBITALS*2;
			for (size_t i=0;i<dof;i++) for (size_t j=i+1;j<dof;j++)
					jmatrix[i+j*dof]=conj(jmatrix[j+i*dof]);

			for (size_t i=0;i<jmatrix.size();i++) jmatrix[i] *= mp_.J;

		}

		void auxCreateYmatrix(std::vector<RealType>& ymatrix,const
				typename DynVarsType::SpinType& dynVars,size_t site) const
		{
			size_t dof = ORBITALS*2;
			ymatrix[28]=mp_.spinOrbitCoupling;
			ymatrix[35]=mp_.spinOrbitCoupling;

			for (size_t i=0;i<dof;i++) for (size_t j=i+1;j<dof;j++)
				ymatrix[i+j*dof]=std::conj(ymatrix[j+i*dof]);
		}


// 		template<typename GreenFunctionType>
// 		RealType calcNumber(GreenFunctionType& greenFunction) const
// 		{
// 			RealType sum=0;
// 			for (size_t i=0;i<hilbertSize_;i++) {
// 				sum += real(greenFunction(i,i));
// 			}
// 			return sum;
// 		}


		
		RealType calcKinetic(const DynVarsType& dynVars,
				      const std::vector<RealType>& eigs) const
		{
			RealType sum = 0;
			//const psimag::Matrix<ComplexType>& matrix = matrix_;
// 			for (size_t lambda=0;lambda<hilbertSize_;lambda++) {
// 				RealType tmp2=0.0;
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

		const EngineParamsType& engineParams_;
		const ParametersModelType& mp_;
		const GeometryType& geometry_;
		DynVarsType dynVars_;
		size_t hilbertSize_;
		AdjustmentsType adjustments_;
		ProgressIndicatorType progress_;
		//RandomNumberGeneratorType& rng_;
		SpinOperationsType spinOperations_;
		ObservablesStoredType observablesStored_;
	}; // PnictidesTwoOrbitals

	template<typename EngineParamsType,
	         typename ParametersModelType,
	         typename GeometryType,
	         typename ConcurrencyType>
	std::ostream& operator<<(std::ostream& os,
	                         const DmsMultiOrbital<EngineParamsType,
	                         ParametersModelType,GeometryType,ConcurrencyType>& model)
	{
		os<<"ModelParameters\n";
		os<<model.mp_;
		return os;
	}
} // namespace Spf

/*@}*/
#endif // DMS_MULIORBITAL_H
