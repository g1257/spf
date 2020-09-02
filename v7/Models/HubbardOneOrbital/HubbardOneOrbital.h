
/** \ingroup SPF */
/*@{*/

/*! \file HubbardOneOrbital
 *
 *  HubbardOneOrbital model
 *
 */
#ifndef HUBBARD_ONE_ORBITAL_H
#define HUBBARD_ONE_ORBITAL_H
#include "HubbardOneOrbitalFields.h"
#include "Random48.h"
#include "ProgressIndicator.h"
#include "../../ClassicalFields/ContVarFiniteOperations.h"
#include "ModelBase.h"
#include "../PnictidesMultiOrbitals/ThreeOrbitalTerms.h"
//#include "HubbardOneOrbitalObsStored.h"
#include "Conductance.h"
#include "CrsMatrix.h"
#include "ParametersHubbardOneOrbital.h"

namespace Spf {
	template<typename EngineParamsType,
	           typename GeometryType>
	class HubbardOneOrbital : public ModelBase<Spin<
	  typename EngineParamsType::RealType>,
	  EngineParamsType,
	  GeometryType> {

		typedef typename EngineParamsType::RealType RealType;
		typedef PsimagLite::CrsMatrix<RealType> SparseMatrixType;
		typedef typename GeometryType::PairType PairType;
		typedef PsimagLite::ProgressIndicator ProgressIndicatorType;

	public:

		typedef std::complex<RealType> ComplexType;
		typedef PsimagLite::Matrix<RealType> MatrixType;
		typedef typename EngineParamsType::IoInType IoInType;
		typedef ParametersHubbardOneOrbital<EngineParamsType,IoInType> ParametersModelType;
		typedef HubbardOneOrbitalFields<RealType,GeometryType> DynVarsType;
		typedef ContVarFiniteOperations<GeometryType,RealType> ContVarFiniteOperationsType;
		typedef typename ContVarFiniteOperationsType::DynVarsType ContVarFiniteType;
		typedef typename ContVarFiniteType::PairRealType PairRealType;
//		typedef HubbardOneOrbitalObsStored<ContVarFiniteOperationsType,RealType,
//				ParametersModelType> HubbardOneOrbitalObsStoredType;

		enum {OLDFIELDS,NEWFIELDS};
		enum {SPIN_UP,SPIN_DOWN};
		enum {CHARGE = DynVarsType::CHARGE, MAG = DynVarsType::MAG};

		HubbardOneOrbital(const EngineParamsType& engineParams,
		                  IoInType& io,
		                  const GeometryType& geometry)
		: engineParams_(engineParams),
		  mp_(io,engineParams),
		  geometry_(geometry),
		  dynVars_(geometry.volume(),engineParams),
		  hilbertSize_(2*geometry.volume()),
		  progress_("HubbardOneOrbital"),
		  chargeOperations_(geometry,engineParams.mcWindow,0,PairRealType(0,2)),
		  magOperations_(geometry,engineParams.mcWindow,1,PairRealType(-1,1))
		{}

		DynVarsType& dynVars() { return dynVars_; }

		SizeType totalFlips() const { return geometry_.volume(); }

		void setOperation(ContVarFiniteOperationsType** op,SizeType i)
		{
			assert(i == 0 || i == 1);
			if (i == 0)
				*op = &chargeOperations_;
			else if (i == 1)
				*op = &magOperations_;
			else
				throw PsimagLite::RuntimeError("HubbardOneOrbital::setOperation()\n");
		}

		SizeType hilbertSize() const { return hilbertSize_; }

		RealType deltaDirect(SizeType i,ContVarFiniteOperationsType& ops,SizeType n) const
		{
			const RealType& coupling = (n==DynVarsType::CHARGE) ? mp_.dampingCharge : mp_.dampingMag;
			return ops.deltaDirect(i,coupling);
		}

		template<typename GreenFunctionType,typename SomePackerType>
		void doMeasurements(GreenFunctionType& greenFunction,SizeType iter,SomePackerType& packer)
		{
			packer.pack("iter=",iter);

			RealType temp=greenFunction.calcNumber();
			packer.pack("Number_Of_Electrons=",temp);

			temp=greenFunction.calcElectronicEnergy();
			packer.pack("EnergyElectronic=",temp);

			RealType temp2 = chargeOperations_.sum2() * mp_.dampingCharge;
			packer.pack("EnergyDampingCharge=",temp2);
			temp += temp2;

			temp2 = magOperations_.sum2() * mp_.dampingMag;
			packer.pack("EnergyDampingMag=",temp2);
			temp += temp2;

			temp2 = chargeOperations_.sum();
			packer.pack("AverageCharge=",temp2);

			temp2 = magOperations_.sum();
			packer.pack("AverageMag=",temp2);

			//! total energy = electronic energy + superexchange + phonon energy
			packer.pack("EnergyTotal=",temp);

			packer.pack("Adjustments: mu=",engineParams_.mu);

//			HubbardOneOrbitalObsStored_(dynVars,greenFunction);
		} // doMeasurements

		void createHamiltonian(MatrixType& matrix,SizeType oldOrNewDynVars)
		{
			DynVarsType newDynVars(chargeOperations_.dynVars2(),magOperations_.dynVars2());

			 if (oldOrNewDynVars==NEWFIELDS) createHamiltonian(newDynVars,matrix);
			 else createHamiltonian(dynVars_,matrix);
		}

		void createHsparse(SparseMatrixType& sparseMatrix,SizeType oldOrNewDynVars)
		{
			// ALL THIS IS VERY INEFFICIENT
			// FIXME, NEEDS TO WRITE THIS FROM SCRATCH!!!!
			MatrixType matrix(hilbertSize_,hilbertSize_);
			createHamiltonian(matrix,oldOrNewDynVars);
			fullMatrixToCrsMatrix(sparseMatrix,matrix);
		}

		void setTpemThings(RealType&, RealType&, PsimagLite::Vector<SizeType>::Type&) const
		{
			throw PsimagLite::RuntimeError("setTpemThings unimplemented\n");
		}

		template<typename SomeOutputType>
		void finalize(SomeOutputType&)
		{
//			HubbardOneOrbitalObsStored_.finalize(fout);
		}

		template<typename EngineParamsType2,
		         typename GeometryType2>
		friend std::ostream& operator<<(std::ostream& os,
		  const HubbardOneOrbital<EngineParamsType2,
		                               GeometryType2>& model);

	private:

		void createHamiltonian(const DynVarsType&,
		                       MatrixType& matrix,
		                       const RealType* = 0) const
		{
			SizeType volume = geometry_.volume();
			SizeType dof = 2; // the 2 comes because of the spin

			for (SizeType gamma1=0;gamma1<matrix.n_row();gamma1++)
				for (SizeType p = 0; p < matrix.n_col(); p++)
					matrix(gamma1,p)=0;

			for (SizeType p = 0; p < volume; p++) {
				for (SizeType spin1=0;spin1<dof;spin1++) {
					matrix(p+spin1*volume,p+spin1*volume) =
					        interaction(p,spin1) + mp_.potentialV[p];
					for (SizeType j = 0; j <  geometry_.z(1); j++) {
						if (j%2!=0) continue;
						PairType tmpPair = geometry_.neighbor(p,j);
						SizeType k = tmpPair.first;

						matrix(p+spin1*volume,k+spin1*volume) = mp_.hopping;
						matrix(k+spin1*volume,p+spin1*volume) =
						        PsimagLite::conj(matrix(p+spin1*volume,k+spin1*volume));

					}
				}
			}
		}

		RealType interaction(SizeType i,SizeType spin) const
		{
			RealType ni = dynVars_.getField(CHARGE).value[i] * mp_.interactionCharge;
			RealType mi = dynVars_.getField(MAG).value[i] * mp_.interactionMag;
			return (spin==SPIN_UP) ? ni+mi : ni-mi;
		}

		const EngineParamsType& engineParams_;
		ParametersModelType mp_;
		const GeometryType& geometry_;
		DynVarsType dynVars_;
		SizeType hilbertSize_;
		ProgressIndicatorType progress_;
		ContVarFiniteOperationsType chargeOperations_;
		ContVarFiniteOperationsType magOperations_;
//		HubbardOneOrbitalObsStoredType HubbardOneOrbitalObsStored_;
	}; // HubbardOneOrbital

	template<typename EngineParamsType,
	         typename GeometryType>
	std::ostream& operator<<(std::ostream& os,
	                         const HubbardOneOrbital<
	                           EngineParamsType,
	                           GeometryType>& model)
	{
		os<<"ModelParameters\n";
		os<<model.mp_;
		return os;
	}
} // namespace Spf

/*@}*/
#endif // HUBBARD_ONE_ORBITAL_H
