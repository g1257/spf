
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
#include "../../ClassicalFields/SpinOperations.h"
#include "ModelBase.h"
#include "ParametersDmsMultiOrbital.h"
#include "DmsMultiOrbitalObsStored.h"
#include "FakeParams.h"
#include "DmsConductance.h"

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
	typedef DmsMultiOrbitalObsStored<SpinOperationsType,
	                                 ComplexType,
	                                 ParametersModelType,
	                                 EngineParamsType> DmsMultiOrbitalObsStoredType;

	enum {OLDFIELDS,NEWFIELDS};

	enum {SPIN_UP, SPIN_DOWN};

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
	{
		if (mp_.potentialV.size() != geometry.volume()) {
			PsimagLite::String str("DmsMultiOrbital: potentialV.size= ");
			str += ttos(mp_.potentialV.size());
			str += " expecting " + ttos(geometry.volume()) + "\n";
			throw PsimagLite::RuntimeError(str);
		}
	}

	DynVarsType& dynVars() { return dynVars_; }

	SizeType totalFlips() const { return geometry_.volume(); }

	void setOperation(SpinOperationsType** op,SizeType i)
	{
		assert(i == 0);
		*op = &spinOperations_;
	}

	SizeType hilbertSize() const { return hilbertSize_; }

	RealType deltaDirect(SizeType i,const SpinOperationsType& ops,int n) const
	{
		assert(n == 0);
		assert(n == 0);
		RealType x = ops.deltaDirect(i,mp_.jafNn,mp_.jafNnn);
		x += ops.deltaMagneticField(i,mp_.magneticField);
		x += ops.deltaDmInteraction(i,mp_.dmNn,mp_.dmNnn);
		return x;
	}

	void set(typename DynVarsType::SpinType& dynVars)
	{
		spinOperations_.set(dynVars);
	}

	template<typename GreenFunctionType,typename SomePackerType>
	void doMeasurements(
	        GreenFunctionType& greenFunction,
	        SizeType iter,
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

		if (mp_.dmNn.size() > 0 || mp_.dmNnn.size() > 0) {
			temp2=spinOperations_.calcDmTerm(dynVars,mp_.dmNn,mp_.dmNnn);
			packer.pack("DmTerm=",temp2);
			temp += temp2;
		}

		// total energy = electronic energy + superexchange + phonon energy
		packer.pack("TotalEnergy=",temp);

		//s="Action=";

		//s="Number_Of_Holes=";

		packer.pack("Adjustments: mu=",engineParams_.mu);

		temp = spinOperations_.calcMag(dynVars,mp_.modulus);
		packer.pack("Mag2=",temp);

		if (engineParams_.options.find("conductance") != PsimagLite::String::npos) {

			DmsConductance<MatrixType, GreenFunctionType> conductance;
			temp = conductance(greenFunction,
			                   engineParams_.mu,
			                   geometry_.dim(),
			                   geometry_.length(),
			                   1000,
			                   1e-6);
			packer.pack("Conductance=",temp);
		}

		DmsMultiOrbitalObsStored_(dynVars,greenFunction);
	} // doMeasurements

	void createHamiltonian(
	        PsimagLite::Matrix<ComplexType>& matrix,
	        SizeType oldOrNewDynVars)
	{
		SpinType* dynVarsPtr = 0;
		dynVars_.getField(&dynVarsPtr,0);
		const SpinType& dynVars = *dynVarsPtr;

		if (oldOrNewDynVars==NEWFIELDS)
			createHamiltonian(spinOperations_.dynVars2(),matrix,mp_.J);
		else
			createHamiltonian(dynVars,matrix,mp_.J);
	}

	void createHsparse(SparseMatrixType& sparseMatrix,SizeType oldOrNewDynVars)
	{
		// ALL THIS IS VERY INEFFICIENT
		// FIXME, NEEDS TO WRITE THIS FROM SCRATCH!!!!
		MatrixType matrix(hilbertSize_,hilbertSize_);
		createHamiltonian(matrix,oldOrNewDynVars);
		fullMatrixToCrsMatrix(sparseMatrix,matrix);
	}

	void setTpemThings(RealType& a,
	                   RealType& b,
	                   PsimagLite::Vector<SizeType>::Type& support) const
	{
		{
			MatrixType matrix(hilbertSize_,hilbertSize_);
			FakeParams fakeParams("none",343313,engineParams_.options);
			SpinType tmpDynVars(geometry_.volume(),fakeParams);
			createHamiltonian(tmpDynVars,matrix,mp_.J);
			typename PsimagLite::Vector<RealType>::Type e(matrix.n_row());
			diag(matrix,e,'N');

			this->setTpemAandB(a,b,e[0],e[e.size()-1]);
		}

		{
			// setup support:
			// to make this work for mp.J==0 we need to set
			RealType J = 0.5;
			MatrixType matrix(hilbertSize_,hilbertSize_);
			FakeParams fakeParams("none",343313,engineParams_.options);
			SpinType tmpDynVars(geometry_.volume(),fakeParams);
			createHamiltonian(tmpDynVars,matrix,J);
			SpinOperationsType spinOps(geometry_,engineParams_);
			spinOps.set(tmpDynVars);
			SizeType site = 0;
			spinOps. makeChange(site,0.2,0.1);

			MatrixType matrix2(hilbertSize_,hilbertSize_);
			createHamiltonian(spinOps.dynVars2(),matrix2,J);

			this->setTpemSupport(support,matrix,matrix2,site);
		}
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
	                       MatrixType& matrix,
	RealType J) const
	{
		SizeType volume = geometry_.volume();
		SizeType norb = mp_.orbitals;
		SizeType dof = norb * 2; // the 2 comes because of the spin
		typename PsimagLite::Vector<ComplexType>::Type jmatrix(dof*dof,0);
		typename PsimagLite::Vector<RealType>::Type ymatrix(dof*dof,0);

		for (SizeType gamma1=0;gamma1<matrix.n_row();gamma1++)
			for (SizeType p = 0; p < matrix.n_col(); p++)
				matrix(gamma1,p)=0;

		for (SizeType p = 0; p < volume; p++) {
			RealType modulus = mp_.modulus[p];

			if (norb == 3) {
				auxCreateJmatrix(jmatrix,dynVars,p,J);
				auxCreateYmatrix(ymatrix,dynVars,p);
			} else {
				auxCreateJmatrix1(jmatrix,dynVars,p,J);
			}

			for (SizeType gamma1=0;gamma1<dof;gamma1++) {
				//! Term B (n_iup - n_idown)
				SizeType spin1 = SizeType(gamma1/norb);
				RealType magField = (spin1==SPIN_UP) ? mp_.magneticField :
				                                       -mp_.magneticField;
				matrix(p+gamma1*volume,p+gamma1*volume) =
				        real(jmatrix[gamma1+dof*gamma1])*modulus +
				        mp_.potentialV[p] + magField +
				        ymatrix[gamma1+dof*gamma1];

				for (SizeType j = 0; j <  geometry_.z(1); j++) {
					//if (j%2!=0) continue;
					PairType tmpPair = geometry_.neighbor(p,j);
					SizeType k = tmpPair.first;
					SizeType dir = static_cast<SizeType>(tmpPair.second/2); //geometry_.scalarDirection(p,k);
					for (SizeType gamma2=0;gamma2<dof;gamma2++) {
						SizeType index = gamma1+gamma2*dof+dir*dof*dof;
						assert(index < mp_.hoppings.size());
						matrix(p+gamma1*volume,k+gamma2*volume) = mp_.hoppings[index];
					}
				}

				for (SizeType j = 0; j <  geometry_.z(2); j++) {
					PairType tmpPair = geometry_.neighbor(p,j,2);
					SizeType k = tmpPair.first;
					SizeType dir = tmpPair.second;
					for (SizeType gamma2=0;gamma2<dof;gamma2++) {
						SizeType index = gamma1+gamma2*dof+dir*dof*dof;
						assert(index < mp_.hoppings.size());
						matrix(p+gamma1*volume,k+gamma2*volume) = mp_.hoppings[index];
					}
				}

				for (SizeType gamma2=0;gamma2<dof;gamma2++) {
					if (gamma1 == gamma2) continue;
					matrix(p+gamma1*volume,p + gamma2*volume) =
					        jmatrix[gamma1+dof*gamma2] * modulus;
				}
			}
		}
	}

	void auxCreateJmatrix1(typename PsimagLite::Vector<ComplexType>::Type& jmatrix,
	                       const typename DynVarsType::SpinType& dynVars,
	                       SizeType site,
	RealType J) const
	{
		jmatrix[0]=cos(dynVars.theta[site]);
		jmatrix[1]=ComplexType(sin(dynVars.theta[site])*cos(dynVars.phi[site]),
		                      -sin(dynVars.theta[site])*sin(dynVars.phi[site]));
		jmatrix[3]= -cos(dynVars.theta[site]);
		jmatrix[2]=conj(jmatrix[1]);

		for (SizeType i=0;i<jmatrix.size();i++) jmatrix[i] *= J;
	}

	void auxCreateJmatrix(typename PsimagLite::Vector<ComplexType>::Type& jmatrix,const
	                      typename DynVarsType::SpinType& dynVars,
	                      SizeType site,
	RealType) const
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

		SizeType dof = mp_.orbitals*2;
		for (SizeType i=0;i<dof;i++) for (SizeType j=i+1;j<dof;j++)
			jmatrix[i+j*dof]=conj(jmatrix[j+i*dof]);

		for (SizeType i=0;i<jmatrix.size();i++) jmatrix[i] *= mp_.J;

	}

	void auxCreateYmatrix(typename PsimagLite::Vector<RealType>::Type& ymatrix,const
	                      typename DynVarsType::SpinType&,SizeType) const
	{
		SizeType dof = mp_.orbitals*2;
		ymatrix[28]=mp_.spinOrbitCoupling;
		ymatrix[35]=mp_.spinOrbitCoupling;

		for (SizeType i=0;i<dof;i++) for (SizeType j=i+1;j<dof;j++)
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
	SizeType hilbertSize_;
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

