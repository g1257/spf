
/** \ingroup SPF */
/*@{*/

/*! \file PhononsTwoOrbitals.h
 *
 *  PhononsTwoOrbitals model (JH infinity)
 *
 */
#ifndef PHONONS_2ORB_H
#define PHONONS_2ORB_H
#include "PhononsTwoOrbitalsFields.h"
#include "Random48.h"
#include "ProgressIndicator.h"
#include "Adjustments.h"
#include "../../ClassicalFields/SpinOperations.h"
#include "../../ClassicalFields/PhononOperations.h"
#include "ModelBase.h"
#include "ParametersPhononsTwoOrbitals.h"

namespace Spf {
	template<typename EngineParamsType,
	         typename GeometryType>
	class PhononsTwoOrbitals : public ModelBase<
	  Spin<typename EngineParamsType::RealType>,
	  EngineParamsType,
	  GeometryType> {

		typedef typename EngineParamsType::RealType RealType;
		typedef std::complex<RealType> ComplexType;
		typedef PsimagLite::CrsMatrix<ComplexType> SparseMatrixType;
		typedef typename GeometryType::PairType PairType;
		typedef PsimagLite::ProgressIndicator ProgressIndicatorType;
		typedef Adjustments<EngineParamsType> AdjustmentsType;
		typedef PhononsTwoOrbitals<EngineParamsType,
		                           GeometryType> ThisType;

		static const SizeType nbands_ = 2;

	public:

		typedef PsimagLite::Matrix<ComplexType> MatrixType;
		typedef typename EngineParamsType::IoInType IoInType;
		typedef ParametersPhononsTwoOrbitals<EngineParamsType,IoInType> ParametersModelType;
		typedef PhononsTwoOrbitalsFields<RealType,GeometryType> DynVarsType;
		typedef typename DynVarsType::SpinType SpinType;
		typedef typename DynVarsType::PhononType PhononType;
		typedef typename DynVarsType::SpinOperationsType SpinOperationsType;
		typedef typename DynVarsType::PhononOperationsType PhononOperationsType;

		enum {OLDFIELDS,NEWFIELDS};

		PhononsTwoOrbitals(EngineParamsType& engineParams,
		                   IoInType& io,
		                   const GeometryType& geometry)
		: engineParams_(engineParams),
		  mp_(io,engineParams),
		  geometry_(geometry),
		  dynVars_(geometry.volume(),engineParams),
		  hilbertSize_(nbands_*geometry_.volume()), // there's no spin here
		  adjustments_(engineParams),
		  progress_("PhononsTwoOrbitals"),
		  spinOperations_(geometry_,engineParams_),
		  phononOperations_(geometry_,
		                    engineParams.mcWindow,
		                    2,
		                    engineParams.options)
		{}

		DynVarsType& dynVars() { return dynVars_; }

		void setOperation(SpinOperationsType** op,SizeType i)
		{
			assert(i == 0);
			*op = &spinOperations_;
		}

		void setOperation(PhononOperationsType** op,SizeType i)
		{
			assert(i == 1);
			*op = &phononOperations_;
		}

		SizeType totalFlips() const { return geometry_.volume(); }

		SizeType hilbertSize() const { return hilbertSize_; }

		RealType deltaDirect(SizeType i,SpinOperationsType& ops,int n) const
		{
			assert(n == 0);
			return ops.deltaDirect(i,mp_.jaf,0);
		}

		RealType deltaDirect(SizeType,PhononOperationsType&,int n) const
		{
			assert(n == 1);
			return 0.0;
		}

		template<typename GreenFunctionType,typename SomePackerType>
		void doMeasurements(GreenFunctionType& greenFunction,SizeType iter,SomePackerType& packer)
		{
			SpinType* dynVarsPtr = 0;
			dynVars_.getField(&dynVarsPtr,0);
			const SpinType& spinPart = *dynVarsPtr;

			packer.pack("iter=",iter);

			RealType temp=calcNumber(greenFunction);
			packer.pack("Number_Of_Electrons=",temp);

			//s = "rankGlobal=";

			temp=greenFunction.calcElectronicEnergy();
			//temp=calcElectronicEnergy(greenFunction);
			packer.pack("Electronic Energy=",temp);

			RealType temp2=spinOperations_.calcSuperExchange(spinPart, mp_.jaf);
			packer.pack("Superexchange=",temp2);

			temp += temp2;

			// total energy = electronic energy + superexchange + phonon energy
			packer.pack("TotalEnergy-FIXME-ADD-PHONON-PART=",temp);
			//s="Action=";

			//s="Number_Of_Holes=";

			packer.pack("Adjustments: mu=",engineParams_.mu);

			temp = spinOperations_.calcMag(spinPart);
			packer.pack("Mag2=",temp);

// 			temp=calcKinetic(dynVars_,eigs);
// 			s ="KineticEnergy="+ttos(temp);
// 			progress_.printline(s,fout);

			//storedObservables_.doThem();
		} // doMeasurements

		void createHamiltonian(PsimagLite::Matrix<ComplexType>& matrix,SizeType oldOrNewDynVars)
		{

			PhononType* phononPartPtr = 0;
			dynVars_.getField(&phononPartPtr,1);
			const PhononType& phononPart = *phononPartPtr;

			DynVarsType newDynVars(spinOperations_.dynVars2(),phononPart);

			 if (oldOrNewDynVars==NEWFIELDS) createHamiltonian(newDynVars,matrix);
			 else createHamiltonian(dynVars_,matrix);
		}

		void createHsparse(SparseMatrixType& sparseMatrix,SizeType oldOrNewDynVars)
		{
			// ALL THIS IS VERY INEFFICIENT
			// FIXME, NEEDS TO WRITE THIS FROM SCRATCH!!!!
			MatrixType matrix;
			createHamiltonian(matrix,oldOrNewDynVars);
			fullMatrixToCrsMatrix(sparseMatrix,matrix);
		}

		void setTpemThings(RealType&, RealType&, PsimagLite::Vector<SizeType>::Type&) const
		{
			throw std::runtime_error("You can't run this model with TPEM yet (sorry)\n");
		}

		void adjustChemPot(const typename PsimagLite::Vector<RealType>::Type& eigs)
		{
			if (engineParams_.carriers==0) return;
			try {
				engineParams_.mu = adjustments_.adjChemPot(eigs);
			} catch (std::exception& e) {
				std::cerr<<e.what()<<"\n";
			}

		}

		template<typename SomeOutputType>
		void finalize(SomeOutputType&)
		{
		}

		template<typename EngineParamsType2,
		         typename GeometryType2>
		friend std::ostream& operator<<(std::ostream& os,
		  const PhononsTwoOrbitals<
		           EngineParamsType2,
		           GeometryType2>& model);

		private:

		void createHamiltonian(DynVarsType&,MatrixType& matrix) const
		{
			SizeType volume = geometry_.volume();

			const SpinType* spinPartPtr = 0;
			dynVars_.getField(&spinPartPtr,0);
			const SpinType& spinPart = *spinPartPtr;

			const PhononType* phononPartPtr = 0;
			dynVars_.getField(&phononPartPtr,1);
			const PhononType& phononPart = *phononPartPtr;

			for (SizeType gamma1=0;gamma1<matrix.n_row();gamma1++)
				for (SizeType p = 0; p < matrix.n_col(); p++)
					matrix(gamma1,p)=0;

			for (SizeType p = 0; p < volume; p++) {
				RealType phonon_q1=phononOperations_.calcPhonon(p,phononPart,0);
				RealType phonon_q2=phononOperations_.calcPhonon(p,phononPart,1);
				RealType phonon_q3=phononOperations_.calcPhonon(p,phononPart,2);
				matrix(p,p) = mp_.phononSpinCoupling[0]*phonon_q1+
						mp_.phononSpinCoupling[2]*phonon_q3+
						mp_.potential[p];
				matrix(p+volume,p+volume) = -mp_.phononSpinCoupling[2]*phonon_q3+
						mp_.phononSpinCoupling[0]*phonon_q1+mp_.potential[p];
				matrix(p,p+volume) = (mp_.phononSpinCoupling[1]*phonon_q2);
				matrix(p+volume,p) = conj(matrix(p,p+volume));

				for (SizeType j = 0; j < geometry_.z(1); j++) {	/* hopping elements, n-n only */
					PairType tmpPair = geometry_.neighbor(p,j);
					SizeType col = tmpPair.first;
					SizeType dir = tmpPair.second; // int(j/2);

					RealType tmp=cos(0.5*spinPart.theta[p])*cos(0.5*spinPart.theta[col]);
					RealType tmp2=sin(0.5*spinPart.theta[p])*sin(0.5*spinPart.theta[col]);
					ComplexType S_ij=ComplexType(tmp+tmp2*cos(spinPart.phi[p]-spinPart.phi[col]),
						-tmp2*sin(spinPart.phi[p]-spinPart.phi[col]));

					matrix(p,col) = -mp_.hoppings[0+0*2+dir*4] * S_ij;
					matrix(col,p) = conj(matrix(p,col));

					matrix(p, col+volume)= -mp_.hoppings[0+1*2+dir*4] * S_ij;
					matrix(col+volume,p)=conj(matrix(p,col+volume));

					matrix(p+volume,col) =  -mp_.hoppings[1+0*2+dir*4] * S_ij;
					matrix(col,p+volume) =  conj(matrix(p+volume,col));

					matrix(p+volume,col+volume) =  -mp_.hoppings[1+1*2+dir*4] * S_ij;
					matrix(col+volume,p+volume) = conj(matrix(p+volume,col+volume));
				}
			}
		}

		template<typename GreenFunctionType>
		RealType calcNumber(GreenFunctionType& greenFunction) const
		{
			RealType sum=0;
			for (SizeType i=0;i<hilbertSize_;i++) {
				sum += real(greenFunction(i,i));
			}
			return sum;
		}

		RealType calcKinetic(const DynVarsType& dynVars,
		                     const typename PsimagLite::Vector<RealType>::Type& eigs) const
		{
			RealType sum = 0;
			//constPsimagLite::Matrix<ComplexType>& matrix = matrix_;
// 			for (SizeType lambda=0;lambda<hilbertSize_;lambda++) {
// 				RealType tmp2=0.0;
// 				for (SizeType i=0;i<geometry_.volume();i++) {
// 					for (SizeType k=0;k<geometry.z(1);k++) {
// 						SizeType j=geometry.neighbor(i,k).first;
// 						for (SizeType spin=0;spin<2;spin++) {
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
		ParametersModelType mp_;
		const GeometryType& geometry_;
		DynVarsType dynVars_;
		SizeType hilbertSize_;
		AdjustmentsType adjustments_;
		ProgressIndicatorType progress_;
		//RandomNumberGeneratorType& rng_;
		SpinOperationsType spinOperations_;
		PhononOperationsType phononOperations_;
	}; // PhononsTwoOrbitals

	template<typename EngineParamsType,
	         typename GeometryType>
	std::ostream& operator<<(std::ostream& os,
      const PhononsTwoOrbitals<
                  EngineParamsType,
                  GeometryType>& model)
	{
		os<<"ModelParameters\n";
		os<<model.mp_;
		return os;
	}
} // namespace Spf

/*@}*/
#endif // PHONONS_2ORB_H
