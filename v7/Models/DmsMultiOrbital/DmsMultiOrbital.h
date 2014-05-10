
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
#include "ParametersDmsMultiOrbital.h"
#include "DmsMultiOrbitalObsStored.h"

namespace Spf {

template<typename EngineParamsType,typename GeometryType>
class DmsMultiOrbital : public ModelBase<Spin<typename EngineParamsType::RealType>,
                                         EngineParamsType,
                                         GeometryType> {

	typedef typename EngineParamsType::RealType RealType;
	typedef std::complex<RealType> ComplexType;
	typedef PsimagLite::CrsMatrix<ComplexType> SparseMatrixType;
	typedef typename GeometryType::PairType PairType;
	typedef PsimagLite::ProgressIndicator ProgressIndicatorType;
	typedef Adjustments<EngineParamsType> AdjustmentsType;
	typedef DmsMultiOrbital<EngineParamsType,GeometryType> ThisType;

public:

	typedef PsimagLite::Matrix<ComplexType> MatrixType;
	typedef typename EngineParamsType::IoInType IoInType;
	typedef ParametersDmsMultiOrbital<EngineParamsType,IoInType> ParametersModelType;
	typedef DmsMultiOrbitalFields<RealType,GeometryType> DynVarsType;
	typedef typename DynVarsType::SpinType SpinType;
	typedef typename DynVarsType::SpinOperationsType SpinOperationsType;
	typedef DmsMultiOrbitalObsStored<SpinOperationsType,ComplexType,
	ParametersModelType,
	EngineParamsType> DmsMultiOrbitalObsStoredType;

	enum {OLDFIELDS,NEWFIELDS};

	DmsMultiOrbital(const EngineParamsType& engineParams,
	                IoInType& io,
	                const GeometryType& geometry)
	    : engineParams_(engineParams),
	      mp_(io,engineParams),
	      geometry_(geometry),
	      dynVars_(geometry.volume(),engineParams,mp_.modulus),
	      hilbertSize_(2*mp_.orbitals*geometry.volume()),
	      adjustments_(engineParams),
	      progress_("PnictidesTwoOrbitals"),
	      spinOperations_(geometry,engineParams),
	      DmsMultiOrbitalObsStored_(spinOperations_,geometry,mp_,engineParams_)
	{}

	DynVarsType& dynVars() { return dynVars_; }

	size_t totalFlips() const { return geometry_.volume(); }

	void setOperation(SpinOperationsType** op,size_t i)
	{
		assert(i == 0);
		*op = &spinOperations_;
	}

	size_t hilbertSize() const { return hilbertSize_; }

	RealType deltaDirect(size_t i,const SpinOperationsType& ops,int n) const
	{
		assert(n == 0);
		return 0.0;
	}

	RealType integrationMeasure(size_t i,const SpinOperationsType& ops,int n)
	{
		assert(n == 0);
		return ops.sineUpdate(i);
	}

	void set(typename DynVarsType::SpinType& dynVars)
	{
		spinOperations_.set(dynVars);
	}

	template<typename GreenFunctionType,typename SomePackerType>
	void doMeasurements(
	        GreenFunctionType& greenFunction,
	        size_t iter,
	        SomePackerType& packer)
	{
		SpinType* dynVarsPtr = 0;
		dynVars_.getField(&dynVarsPtr,0);
		const SpinType& dynVars = *dynVarsPtr;

		packer.pack("iter=",iter);

		RealType temp=greenFunction.calcNumber();
		packer.pack("Number_Of_Electrons=",temp);

		temp=greenFunction.calcElectronicEnergy();
		packer.pack("Electronic Energy=",temp);

		RealType temp2=spinOperations_.calcSuperExchange(dynVars,mp_.jafNn);
		packer.pack("Superexchange=",temp2);

		temp += temp2;
		if (mp_.jafNnn!=0) {
			temp2=spinOperations_.directExchange2(dynVars,mp_.jafNnn);
			packer.pack("Superexchange2=",temp2);
			temp += temp2;
		}

		// total energy = electronic energy + superexchange + phonon energy
		packer.pack("TotalEnergy=",temp);

		//s="Action=";

		//s="Number_Of_Holes=";

		packer.pack("Adjustments: mu=",engineParams_.mu);

		temp = spinOperations_.calcMag(dynVars,mp_.modulus);
		packer.pack("Mag2=",temp);


		// 			temp=calcKinetic(dynVars_,eigs);
		// 			s ="KineticEnergy="+ttos(temp);
		// 			progress_.printline(s,fout);

		DmsMultiOrbitalObsStored_(dynVars,greenFunction);
	} // doMeasurements

	void createHamiltonian(
	        PsimagLite::Matrix<ComplexType>& matrix,
	        size_t oldOrNewDynVars)
	{
		SpinType* dynVarsPtr = 0;
		dynVars_.getField(&dynVarsPtr,0);
		const SpinType& dynVars = *dynVarsPtr;

		if (oldOrNewDynVars==NEWFIELDS)
			createHamiltonian(spinOperations_.dynVars2(),matrix);
		else
			createHamiltonian(dynVars,matrix);
	}

	void createHsparse(SparseMatrixType& sparseMatrix,size_t oldOrNewDynVars)
	{
		// ALL THIS IS VERY INEFFICIENT
		// FIXME, NEEDS TO WRITE THIS FROM SCRATCH!!!!
		MatrixType matrix(hilbertSize_,hilbertSize_);
		createHamiltonian(matrix,oldOrNewDynVars);
		fullMatrixToCrsMatrix(sparseMatrix,matrix);
	}

	void setTpemThings(RealType& a,
	                   RealType& b,
	                   PsimagLite::Vector<size_t>::Type& support) const
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
	void finalize(SomeOutputType& fout)
	{
		DmsMultiOrbitalObsStored_.finalize(fout);
	}

	template<typename EngineParamsType2,typename GeometryType2>
	friend std::ostream& operator<<(std::ostream& os,
	                                const DmsMultiOrbital<EngineParamsType2,
	                                GeometryType2>& model);

private:

	void createHamiltonian(const typename DynVarsType::SpinType& dynVars,
	                       MatrixType& matrix) const
	{
		size_t volume = geometry_.volume();
		size_t norb = mp_.orbitals;
		size_t dof = norb * 2; // the 2 comes because of the spin
		typename PsimagLite::Vector<ComplexType>::Type jmatrix(dof*dof,0);
		typename PsimagLite::Vector<RealType>::Type ymatrix(dof*dof,0);

		for (size_t gamma1=0;gamma1<matrix.n_row();gamma1++)
			for (size_t p = 0; p < matrix.n_col(); p++)
				matrix(gamma1,p)=0;

		for (size_t p = 0; p < volume; p++) {
			RealType modulus = mp_.modulus[p];

			if (norb == 3) {
				auxCreateJmatrix(jmatrix,dynVars,p);
				auxCreateYmatrix(ymatrix,dynVars,p);
			} else {
				auxCreateJmatrix1(jmatrix,dynVars,p);
			}

			for (size_t gamma1=0;gamma1<dof;gamma1++) {

				matrix(p+gamma1*volume,p+gamma1*volume) =
				        real(jmatrix[gamma1+dof*gamma1])*modulus +
				        mp_.potentialV[p] +
				        ymatrix[gamma1+dof*gamma1];

				for (size_t j = 0; j <  geometry_.z(1); j++) {
					//if (j%2!=0) continue;
					PairType tmpPair = geometry_.neighbor(p,j);
					size_t k = tmpPair.first;
					size_t dir = tmpPair.second; //geometry_.scalarDirection(p,k);
					for (size_t gamma2=0;gamma2<dof;gamma2++) {
						SizeType index = gamma1+gamma2*dof+dir*dof*dof;
						assert(index < mp_.hoppings.size());
						matrix(p+gamma1*volume,k+gamma2*volume) = mp_.hoppings[index];
					}
				}

				for (size_t gamma2=0;gamma2<dof;gamma2++) {
					if (gamma1 == gamma2) continue;
					matrix(p+gamma1*volume,p + gamma2*volume) =
					        jmatrix[gamma1+dof*gamma2] * modulus;
				}
			}
		}
	}

	void auxCreateJmatrix1(typename PsimagLite::Vector<ComplexType>::Type& jmatrix,
	                       const typename DynVarsType::SpinType& dynVars,
	                       size_t site) const
	{
		jmatrix[0]=cos(dynVars.theta[site]);
		jmatrix[1]=ComplexType(sin(dynVars.theta[site])*cos(dynVars.phi[site]),
		                      -sin(dynVars.theta[site])*sin(dynVars.phi[site]));
		jmatrix[3]= -cos(dynVars.theta[site]);
		jmatrix[2]=conj(jmatrix[1]);
	}

	void auxCreateJmatrix(typename PsimagLite::Vector<ComplexType>::Type& jmatrix,const
	                      typename DynVarsType::SpinType& dynVars,
	                      size_t site) const
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

		size_t dof = mp_.orbitals*2;
		for (size_t i=0;i<dof;i++) for (size_t j=i+1;j<dof;j++)
			jmatrix[i+j*dof]=conj(jmatrix[j+i*dof]);

		for (size_t i=0;i<jmatrix.size();i++) jmatrix[i] *= mp_.J;

	}

	void auxCreateYmatrix(typename PsimagLite::Vector<RealType>::Type& ymatrix,const
	                      typename DynVarsType::SpinType& dynVars,size_t site) const
	{
		size_t dof = mp_.orbitals*2;
		ymatrix[28]=mp_.spinOrbitCoupling;
		ymatrix[35]=mp_.spinOrbitCoupling;

		for (size_t i=0;i<dof;i++) for (size_t j=i+1;j<dof;j++)
			ymatrix[i+j*dof]=std::conj(ymatrix[j+i*dof]);
	}

	RealType calcKinetic(const DynVarsType& dynVars,
	                     const typename PsimagLite::Vector<RealType>::Type& eigs) const
	{
		return 0.0;
	}

	const EngineParamsType& engineParams_;
	ParametersModelType mp_;
	const GeometryType& geometry_;
	DynVarsType dynVars_;
	size_t hilbertSize_;
	AdjustmentsType adjustments_;
	ProgressIndicatorType progress_;
	SpinOperationsType spinOperations_;
	DmsMultiOrbitalObsStoredType DmsMultiOrbitalObsStored_;
}; // PnictidesTwoOrbitals

template<typename EngineParamsType,
         typename GeometryType>
std::ostream& operator<<(std::ostream& os,
                         const DmsMultiOrbital<EngineParamsType,GeometryType>& model)
{
	os<<"ModelParameters\n";
	os<<model.mp_;
	return os;
}
} // namespace Spf

/*@}*/
#endif // DMS_MULIORBITAL_H
