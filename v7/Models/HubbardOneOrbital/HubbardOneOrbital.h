
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
#include "ContVarFiniteOperations.h"
#include "ModelBase.h"
#include "ThreeOrbitalTerms.h"
//#include "HubbardOneOrbitalObsStored.h"
#include "Conductance.h"
#include "CrsMatrix.h"
#include "ParametersHubbardOneOrbital.h"

namespace Spf {
	template<typename EngineParamsType,
	           typename GeometryType,
	           typename ConcurrencyType_>
	class HubbardOneOrbital : public ModelBase<Spin<
	  typename EngineParamsType::RealType>,
	  EngineParamsType,
	  GeometryType,
	  ConcurrencyType_> {

		typedef typename EngineParamsType::RealType RealType;
		typedef PsimagLite::CrsMatrix<RealType> SparseMatrixType;
		typedef typename GeometryType::PairType PairType;
		typedef PsimagLite::ProgressIndicator ProgressIndicatorType;

	public:

		typedef ConcurrencyType_ ConcurrencyType;
		typedef typename ConcurrencyType::CommType CommType;
		typedef std::complex<RealType> ComplexType;
		typedef PsimagLite::Matrix<RealType> MatrixType;
		typedef typename EngineParamsType::IoInType IoInType;
		typedef ParametersHubbardOneOrbital<EngineParamsType,IoInType> ParametersModelType;
		typedef HubbardOneOrbitalFields<RealType,GeometryType> DynVarsType;
		typedef ContVarFiniteOperations<GeometryType,RealType> ContVarFiniteOperationsType;
		typedef typename ContVarFiniteOperationsType::DynVarsType ContVarFiniteType;
		typedef typename ContVarFiniteType::PairRealType PairRealType;
//		typedef HubbardOneOrbitalObsStored<ContVarFiniteOperationsType,RealType,
//				ParametersModelType,ConcurrencyType> HubbardOneOrbitalObsStoredType;
		
		enum {OLDFIELDS,NEWFIELDS};
		enum {SPIN_UP,SPIN_DOWN};
		
		HubbardOneOrbital(const EngineParamsType& engineParams,
		                  IoInType& io,
		                  const GeometryType& geometry,
		                  ConcurrencyType& concurrency)
		: engineParams_(engineParams),
		  mp_(io,engineParams),
		  geometry_(geometry),
		  concurrency_(concurrency),
		  dynVars_(geometry.volume(),engineParams),
		  hilbertSize_(2*geometry.volume()),
		  progress_("HubbardOneOrbital",0),
		  chargeOperations_(geometry,engineParams.mcWindow.find("Charge")->second,PairRealType(0,2)),
		  magOperations_(geometry,engineParams.mcWindow.find("Mag")->second,PairRealType(-1,1))
//		  HubbardOneOrbitalObsStored_(chargeOperations_,geometry,mp_,2*mp_.numberOfOrbitals,concurrency)
		{
			ProgramGlobals::checkMcWindow(engineParams.mcWindow,"Mag");
			ProgramGlobals::checkMcWindow(engineParams.mcWindow,"Charge");
		}
		
		DynVarsType& dynVars() { return dynVars_; }
		
		size_t totalFlips() const { return geometry_.volume(); }

		void setOperation(ContVarFiniteOperationsType** op,size_t i)
		{
			assert(i == 0 || i == 1);
			if (i == 0)
				*op = &chargeOperations_;
			else if (i == 1)
				*op = &magOperations_;
			else
				throw PsimagLite::RuntimeError("HubbardOneOrbital::setOperation()\n");
		}
		
		size_t hilbertSize() const { return hilbertSize_; }
		
		ConcurrencyType& concurrency() { return concurrency_; }

		RealType deltaDirect(size_t i,ContVarFiniteOperationsType& ops,size_t n) const
		{
			return 0.0;
		}
		
		template<typename GreenFunctionType,typename SomePackerType>
		void doMeasurements(GreenFunctionType& greenFunction,size_t iter,SomePackerType& packer)
		{
			ContVarFiniteType* dynVarsPtr = 0;
			dynVars_.getField(&dynVarsPtr,0);
			const ContVarFiniteType& dynVars = *dynVarsPtr;
			
			packer.pack("iter=",iter);

			RealType temp=greenFunction.calcNumber();
// 			s ="Number_Of_Electrons="+ttos(temp);
			packer.pack("Number_Of_Electrons=",temp);
			
			//s = "rankGlobal=";
			
			temp=greenFunction.calcElectronicEnergy();
// 			s="Electronic Energy="+ttos(temp);
			packer.pack("Electronic Energy=",temp);
			
			RealType temp2=0; //chargeOperations_.calcSuperExchange(dynVars,mp_.jafNn);
			packer.pack("Superexchange=",temp2);
			
			temp += temp2;

			//! total energy = electronic energy + superexchange + phonon energy
			packer.pack("TotalEnergy=",temp);

			packer.pack("Adjustments: mu=",engineParams_.mu);
	
			typename PsimagLite::Vector<RealType>::Type magVector(3,0);
//			chargeOperations_.calcMagVector(magVector,dynVars);
			packer.pack("ClassicalMagnetizationSquared=",magVector*magVector);

			typename PsimagLite::Vector<ComplexType>::Type electronSpinVector(3,0);
			greenFunction.electronSpin(electronSpinVector,1,dynVars.size);
			typename PsimagLite::Vector<ComplexType>::Type combinedVector(3,0);
			combinedVector =  electronSpinVector + magVector;
// 			s="CombinedMagnetizationSquared="+ttos(std::real(combinedVector*combinedVector));
// 			progress_.printline(s,fout);
			packer.pack("CombinedMagnetizationSquared=",
						std::real(combinedVector*combinedVector));

			for (size_t i = 0;i<combinedVector.size();i++) {
// 				s="CombinedMagnetization"+ttos(i)+"="+ttos(combinedVector[i]);
				packer.pack("CombinedMagnetization"+ttos(i)+"=",combinedVector[i]);
// 				progress_.printline(s,fout);
			}

//			HubbardOneOrbitalObsStored_(dynVars,greenFunction);
		} // doMeasurements

		void createHamiltonian(MatrixType& matrix,size_t oldOrNewDynVars)
		{
			DynVarsType newDynVars(chargeOperations_.dynVars2(),magOperations_.dynVars2());

			 if (oldOrNewDynVars==NEWFIELDS) createHamiltonian(newDynVars,matrix);
			 else createHamiltonian(dynVars_,matrix);
		}

		void createHsparse(SparseMatrixType& sparseMatrix,size_t oldOrNewDynVars)
		{
			// ALL THIS IS VERY INEFFICIENT
			// FIXME, NEEDS TO WRITE THIS FROM SCRATCH!!!!
			MatrixType matrix(hilbertSize_,hilbertSize_);
			createHamiltonian(matrix,oldOrNewDynVars);
			fullMatrixToCrsMatrix(sparseMatrix,matrix); 
		}

		struct FakeParams {
			FakeParams(PsimagLite::String dynvarsfile1,int long long randomSeed1)
			: dynvarsfile(dynvarsfile1),randomSeed(randomSeed1) 
			{}

			PsimagLite::String dynvarsfile;
			int long long randomSeed;
		};

		void setTpemThings(RealType& a, RealType& b, PsimagLite::Vector<size_t>::Type& support) const
		{
			throw PsimagLite::RuntimeError("setTpemThings unimplemented\n");
		}

		template<typename SomeOperationsType>
		RealType integrationMeasure(size_t i,SomeOperationsType& ops,int n)
		{
			return 1.0;
		}

		template<typename SomeOutputType>
		void finalize(SomeOutputType& fout,CommType comm)
		{
//			HubbardOneOrbitalObsStored_.finalize(fout,comm);
		}

		template<typename EngineParamsType2,
		         typename GeometryType2,
		         typename ConcurrencyType2>
		friend std::ostream& operator<<(std::ostream& os,
		  const HubbardOneOrbital<EngineParamsType2,
		                               GeometryType2,
		                               ConcurrencyType2>& model);

	private:

		void createHamiltonian(const DynVarsType& dynVars,
		                       MatrixType& matrix,
		                       const RealType* J = 0) const
		{
			size_t volume = geometry_.volume();
			size_t dof = 2; // the 2 comes because of the spin
			typename PsimagLite::Vector<RealType>::Type jmatrix(2*2,0);

			for (size_t gamma1=0;gamma1<matrix.n_row();gamma1++)
				for (size_t p = 0; p < matrix.n_col(); p++)
					matrix(gamma1,p)=0;

			for (size_t p = 0; p < volume; p++) {
				for (size_t spin1=0;spin1<dof;spin1++) {
					matrix(p+spin1*volume,p+spin1*volume) =
					        std::real(jmatrix[spin1+2*spin1]) + mp_.potentialV[p];
					for (size_t j = 0; j <  geometry_.z(1); j++) {
						if (j%2!=0) continue;
						PairType tmpPair = geometry_.neighbor(p,j);
						size_t k = tmpPair.first;

						matrix(p+spin1*volume,k+spin1*volume) = mp_.hopping;
						matrix(k+spin1*volume,p+spin1*volume) = std::conj(matrix(p+spin1*volume,k+spin1*volume));

					}

					for (size_t spin2=0;spin2<2;spin2++) {
						if (spin1==spin2) continue; // diagonal term already taken into account
						matrix(p+spin1*volume,p + spin2*volume)=jmatrix[spin1+2*spin2];
						matrix(p + spin2*volume,p+spin1*volume) =
						        std::conj(matrix(p+spin1*volume,p + spin2*volume));
					}
				}
			}
		}

		const EngineParamsType& engineParams_;
		ParametersModelType mp_;
		const GeometryType& geometry_;
		ConcurrencyType& concurrency_;
		DynVarsType dynVars_;
		size_t hilbertSize_;
		ProgressIndicatorType progress_;
		ContVarFiniteOperationsType chargeOperations_;
		ContVarFiniteOperationsType magOperations_;
//		HubbardOneOrbitalObsStoredType HubbardOneOrbitalObsStored_;
	}; // HubbardOneOrbital

	template<typename EngineParamsType,
	         typename GeometryType,
	         typename ConcurrencyType>
	std::ostream& operator<<(std::ostream& os,
	                         const HubbardOneOrbital<
	                           EngineParamsType,
	                           GeometryType,
	                           ConcurrencyType>& model)
	{
		os<<"ModelParameters\n";
		os<<model.mp_;
		return os;
	}
} // namespace Spf

/*@}*/
#endif // HUBBARD_ONE_ORBITAL_H
