
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
#include "SpinOperations.h"
#include "ModelBase.h"
#include "ThreeOrbitalTerms.h"
#include "PnictidesMultiOrbitalsObsStored.h"
#include "Conductance.h"
#include "CrsMatrix.h"
#include "ParametersPnictidesThreeOrbitals.h"
#include "FakeParams.h"

namespace Spf {
	template<typename EngineParamsType,
	           typename GeometryType>
	class PnictidesMultiOrbitals : public ModelBase<Spin<
	  typename EngineParamsType::RealType>,
	  EngineParamsType,
	  GeometryType> {

		typedef typename EngineParamsType::RealType RealType;
		typedef std::complex<RealType> ComplexType;
		typedef PsimagLite::CrsMatrix<ComplexType> SparseMatrixType;
		//typedef RandomNumberGenerator<RealType> RandomNumberGeneratorType;
		typedef typename GeometryType::PairType PairType;
		typedef PsimagLite::ProgressIndicator ProgressIndicatorType;

	public:

		typedef PsimagLite::Matrix<ComplexType> MatrixType;
		typedef typename EngineParamsType::IoInType IoInType;
		typedef ParametersPnictidesThreeOrbitals<EngineParamsType,IoInType> ParametersModelType;
		typedef PnictidesTwoOrbitalsFields<RealType,GeometryType> DynVarsType;
		typedef typename DynVarsType::SpinType SpinType;
		typedef typename DynVarsType::SpinOperationsType SpinOperationsType;
		typedef ThreeOrbitalTerms<MatrixType,ParametersModelType,
				GeometryType> ThreeOrbitalTermsType;
		typedef PnictidesMultiOrbitalsObsStored<SpinOperationsType,ComplexType,
				ParametersModelType> PnictidesMultiOrbitalsObsStoredType;
		
		enum {OLDFIELDS,NEWFIELDS};
		enum {SPIN_UP,SPIN_DOWN};
		
		PnictidesMultiOrbitals(const EngineParamsType& engineParams,
		                       IoInType& io,
		                       const GeometryType& geometry)
		: engineParams_(engineParams),
		  mp_(io,engineParams),
		  geometry_(geometry),
		  dynVars_(geometry.volume(),engineParams),
		  hilbertSize_(2*mp_.numberOfOrbitals*geometry.volume()),
		  progress_("PnictidesTwoOrbitals"),
		  spinOperations_(geometry,engineParams),
		  threeOrbitalTerms_(mp_,geometry),
		  PnictidesMultiOrbitalsObsStored_(spinOperations_,geometry,mp_,2*mp_.numberOfOrbitals)
		{
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
			RealType x = ops.deltaDirect(i,mp_.jafNn,mp_.jafNnn);
			x += ops.deltaMagneticField(i,mp_.magneticField);
			return x;
		}

		template<typename SomeOperationsType>
		RealType integrationMeasure(SizeType i,SomeOperationsType& ops,int n)
		{
			assert(n == 0);
			return ops.sineUpdate(i);
		}
		
		void set(typename DynVarsType::SpinType& dynVars) { spinOperations_.set(dynVars); }
		
		template<typename GreenFunctionType,typename SomePackerType>
		void doMeasurements(GreenFunctionType& greenFunction,SizeType iter,SomePackerType& packer)
		{
			SpinType* dynVarsPtr = 0;
			dynVars_.getField(&dynVarsPtr,0);
			const SpinType& dynVars = *dynVarsPtr;
			
			packer.pack("iter=",iter);

			RealType temp=greenFunction.calcNumber();
// 			s ="Number_Of_Electrons="+ttos(temp);
			packer.pack("Number_Of_Electrons=",temp);
			
			//s = "rankGlobal=";
			
			temp=greenFunction.calcElectronicEnergy();
// 			s="Electronic Energy="+ttos(temp);
			packer.pack("Electronic Energy=",temp);
			
			RealType temp2=spinOperations_.calcSuperExchange(dynVars,mp_.jafNn);
// 			s="Superexchange="+ttos(temp2);
			packer.pack("Superexchange=",temp2);
// 			progress_.printline(s,fout);
			
			temp += temp2;
			if (mp_.jafNnn!=0) {
				temp2=spinOperations_.directExchange2(dynVars,mp_.jafNnn);
// 				s="Superexchange2="+ttos(temp2);
				packer.pack("Superexchange2=",temp2);
// 				progress_.printline(s,fout);
				temp += temp2;
			}

			//! total energy = electronic energy + superexchange + phonon energy
// 			s="TotalEnergy="+ttos(temp);
			packer.pack("TotalEnergy=",temp);
// 			progress_.printline(s,fout);

			packer.pack("Adjustments: mu=",engineParams_.mu);
	
			typename PsimagLite::Vector<RealType>::Type magVector(3,0);
			spinOperations_.calcMagVector(magVector,dynVars);
// 			s="ClassicalMagnetizationSquared="+ttos(magVector*magVector);
			packer.pack("ClassicalMagnetizationSquared=",magVector*magVector);
// 			progress_.printline(s,fout);

			typename PsimagLite::Vector<ComplexType>::Type electronSpinVector(3,0);
			greenFunction.electronSpin(electronSpinVector,mp_.numberOfOrbitals,dynVars.size);
			typename PsimagLite::Vector<ComplexType>::Type combinedVector(3,0);
			combinedVector =  electronSpinVector + magVector;
// 			s="CombinedMagnetizationSquared="+ttos(std::real(combinedVector*combinedVector));
// 			progress_.printline(s,fout);
			packer.pack("CombinedMagnetizationSquared=",
						std::real(combinedVector*combinedVector));

			for (SizeType i = 0;i<combinedVector.size();i++) {
// 				s="CombinedMagnetization"+ttos(i)+"="+ttos(combinedVector[i]);
				packer.pack("CombinedMagnetization"+ttos(i)+"=",combinedVector[i]);
// 				progress_.printline(s,fout);
			}

			if (engineParams_.options.find("conductance")!=PsimagLite::String::npos) {
				//greenFunction.printMatrix(OLDFIELDS);
				PsimagLite::Matrix<RealType> v
					(greenFunction.hilbertSize(),greenFunction.hilbertSize());
				calcVelocitySquared(greenFunction,v,ProgramGlobals::DIRX);
				typedef Conductance<EngineParamsType,GreenFunctionType> ConductanceType;
				ConductanceType conductance(engineParams_,greenFunction);
// 				s = "ConductanceX=" + ttos(conductance(v));
				packer.pack("ConductanceX=" ,conductance(v));
// 				progress_.printline(s,fout);
				//PsimagLite::Matrix<RealType> vv = v;
				calcVelocitySquared(greenFunction,v,ProgramGlobals::DIRY);
// 				s = "ConductanceY=" + ttos(conductance(v));
				packer.pack("ConductanceY=" ,conductance(v));
// 				progress_.printline(s,fout);
			}

			PnictidesMultiOrbitalsObsStored_(dynVars,greenFunction);
		} // doMeasurements

		void createHamiltonian(MatrixType& matrix,SizeType oldOrNewDynVars)
		{
			SpinType* dynVarsPtr = 0;
			dynVars_.getField(&dynVarsPtr,0);
			const SpinType& dynVars = *dynVarsPtr;

			if (oldOrNewDynVars==NEWFIELDS) {
				createHamiltonian(spinOperations_.dynVars2(),matrix);
//				std::cout<<"--------------\n";
//				std::cout<<spinOperations_.dynVars2();
			} else {
				createHamiltonian(dynVars,matrix);
			}
		}

		void createHsparse(SparseMatrixType& sparseMatrix,SizeType oldOrNewDynVars)
		{
			// ALL THIS IS VERY INEFFICIENT
			// FIXME, NEEDS TO WRITE THIS FROM SCRATCH!!!!
			MatrixType matrix(hilbertSize_,hilbertSize_);
			createHamiltonian(matrix,oldOrNewDynVars);
			fullMatrixToCrsMatrix(sparseMatrix,matrix); 
		}

		void setTpemThings(RealType& a, RealType& b, PsimagLite::Vector<SizeType>::Type& support) const
		{
			{
				MatrixType matrix(hilbertSize_,hilbertSize_);
				FakeParams fakeParams("none",343313,engineParams_.options);
				SpinType tmpDynVars(geometry_.volume(),fakeParams); 
				createHamiltonian(tmpDynVars,matrix);
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
				createHamiltonian(tmpDynVars,matrix,&J);
				SpinOperationsType spinOps(geometry_,engineParams_);
				spinOps.set(tmpDynVars);
				SizeType site = 0;
// 				PsimagLite::Random48<RealType> rng(3493891);
				spinOps. makeChange(site,0.2,0.1);

				MatrixType matrix2(hilbertSize_,hilbertSize_);
				createHamiltonian(spinOps.dynVars2(),matrix2,&J);

				this->setTpemSupport(support,matrix,matrix2,site);
			}
		}

		template<typename SomeOutputType>
		void finalize(SomeOutputType& fout)
		{
			PnictidesMultiOrbitalsObsStored_.finalize(fout);
		}

		template<typename EngineParamsType2,
		         typename GeometryType2>
		friend std::ostream& operator<<(std::ostream& os,
		  const PnictidesMultiOrbitals<EngineParamsType2,
		                               GeometryType2>& model);

	private:

		void createHamiltonian(const typename DynVarsType::SpinType& dynVars,
		                       MatrixType& matrix,
		                       const RealType* J = 0) const
		{
			SizeType volume = geometry_.volume();
			SizeType norb = mp_.numberOfOrbitals;
			SizeType dof = norb * 2; // the 2 comes because of the spin
			typename PsimagLite::Vector<ComplexType>::Type jmatrix(2*2,0);

			for (SizeType gamma1=0;gamma1<matrix.n_row();gamma1++) 
				for (SizeType p = 0; p < matrix.n_col(); p++) 
					matrix(gamma1,p)=0;
			RealType jvalue = (J==0) ? mp_.J : *J;
			for (SizeType p = 0; p < volume; p++) {
				auxCreateJmatrix(jmatrix,dynVars,p,jvalue);
				for (SizeType gamma1=0;gamma1<dof;gamma1++) {
					SizeType spin1 = SizeType(gamma1/norb);
					SizeType orb1 = gamma1 % norb;
					//! Term B (n_iup - n_idown)
					RealType magField = (spin1==SPIN_UP) ? mp_.magneticField : -mp_.magneticField;
					matrix(p+gamma1*volume,p+gamma1*volume) =
						real(jmatrix[spin1+2*spin1]) + mp_.potentialV[p] + magField;
					for (SizeType j = 0; j <  geometry_.z(1); j++) {	
						if (j%2!=0) continue;	
						PairType tmpPair = geometry_.neighbor(p,j);
						SizeType k = tmpPair.first;
						SizeType dir = tmpPair.second; // int(j/2);
						for (SizeType orb2=0;orb2<norb;orb2++) {
							SizeType gamma2 = orb2+spin1*norb; // spin2 == spin1 here
							matrix(p+gamma1*volume,k+gamma2*volume)=
								mp_.hoppings[orb1+orb2*norb+norb*norb*dir];
							matrix(k+gamma2*volume,p+gamma1*volume) = conj(matrix(p+gamma1*volume,k+gamma2*volume));
						}
					}
					//if (geometry.z(p,2)!=4 || geometry.z(p)!=4) throw PsimagLite::RuntimeError("neighbours wrong\n");
					for (SizeType j = 0; j <  geometry_.z(2); j++) {
						if (j%2!=0) continue;	
						PairType tmpPair = geometry_.neighbor(p,j,2);
						SizeType k = tmpPair.first;
						SizeType dir = tmpPair.second;
						//std::cerr<<"Neigbors "<<p<<" "<<k<<"\n";
						for (SizeType orb2=0;orb2<norb;orb2++) {
							SizeType gamma2 = orb2+spin1*norb; // spin2 == spin1 here
							matrix(p+gamma1*volume,k+gamma2*volume)=
								mp_.hoppings[orb1+orb2*norb+norb*norb*dir];
							matrix(k+gamma2*volume,p+gamma1*volume) = conj(matrix(p+gamma1*volume,k+gamma2*volume));
						}
					}
					
					for (SizeType spin2=0;spin2<2;spin2++) {
						if (spin1==spin2) continue; // diagonal term already taken into account
						SizeType gamma2 = orb1+spin2*norb; // orb2 == orb1 here
						matrix(p+gamma1*volume,p + gamma2*volume)=jmatrix[spin1+2*spin2];
						matrix(p + gamma2*volume,p+gamma1*volume) =
								conj(matrix(p+gamma1*volume,p + gamma2*volume));
					}
				}
			}
			threeOrbitalTerms_(matrix);
		}

		void auxCreateJmatrix(typename PsimagLite::Vector<ComplexType>::Type& jmatrix,const
		                      typename DynVarsType::SpinType& dynVars,SizeType site,
		                      const RealType& J) const
		{
			if (!mp_.modulus[site]) {
				for (SizeType i=0;i<jmatrix.size();i++) jmatrix[i] = 0;
				return;
			}

			jmatrix[0]=cos(dynVars.theta[site]);

			jmatrix[1]=ComplexType(sin(dynVars.theta[site])*cos(dynVars.phi[site]),
				sin(dynVars.theta[site])*sin(dynVars.phi[site]));
		
			jmatrix[2]=conj(jmatrix[1]);

			jmatrix[3]= -cos(dynVars.theta[site]);

			for (SizeType i=0;i<jmatrix.size();i++) jmatrix[i] *= J;
		}

		SizeType getSiteAtLayer(SizeType xOrY,SizeType layer,SizeType dir) const
		{
			return (dir==0) ? geometry_.coorToIndex(xOrY,layer) : geometry_.coorToIndex(layer,xOrY);
		}

		template<typename GreenFunctionType>
		void calcVelocitySquared(const GreenFunctionType& gf,
		                     PsimagLite::Matrix<RealType>& v,
		                     SizeType dir) const
		{
			PsimagLite::Matrix<ComplexType> w(v.n_row(),v.n_col());
			SizeType l = geometry_.length();
			for (SizeType layer=0;layer<l;layer++) {
				calcVelocitySquared(gf,w,dir,layer);
			}

			for (SizeType a=0;a<v.n_row();a++) {
				for (SizeType b=0;b<v.n_col();b++) {
					v(a,b) = std::real(std::conj(w(a,b))*w(a,b));
					//if (fabs(v(a,b))>1e-12 && dir==1) std::cout<<"   "<<(a+1)<<"  "<<(b+1)<<"  "<<v(a,b)<<"\n";
				}
			}
			//std::cout<<"----------------------\n";
		}

		template<typename GreenFunctionType>
		void calcVelocitySquared(const GreenFunctionType& gf,
								 PsimagLite::Matrix<ComplexType>& w,
		                         SizeType dir,
								 SizeType startingLayer) const
		{
			SizeType ly = geometry_.length();
			SizeType norb = mp_.numberOfOrbitals;
			SizeType volume = geometry_.volume();

			for (SizeType y=0;y<ly;y++) {
				for (SizeType spin=0;spin<2;spin++) {
					for (SizeType orb1=0;orb1<norb;orb1++) {
						SizeType gamma1 = orb1 + spin*norb;
						SizeType isite = getSiteAtLayer(y,startingLayer,dir);
						SizeType i = isite + gamma1*volume; // starting layer

						for (SizeType y2=0;y2<ly;y2++) {
							for (SizeType orb2=0;orb2<norb;orb2++) {
								SizeType gamma2 = orb2 + spin*norb;
								SizeType endingLayer = startingLayer+1;
								if (endingLayer==ly) endingLayer = 0;
								SizeType jsite = getSiteAtLayer(y2,endingLayer,dir);
								SizeType j = jsite + gamma2*volume; // ending layer

								int dir2 = geometry_.getDirection(isite,jsite);
								if (dir2<0) continue;

								SizeType h = orb1+orb2*norb+norb*norb*dir2;
								RealType hopping = mp_.hoppings[h] 
								    + threeOrbitalTerms_.hopping(isite,dir2,orb1,orb2);
								if (fabs(hopping)<1e-8) continue;
								//std::cerr<<"dir="<<dir<<" isite="<<isite<<" jsite="<<jsite<<" dir2="<<dir2<<" h="<<hopping<<"\n";
								for (SizeType a=0;a<w.n_row();a++) {
									for (SizeType b=0;b<w.n_col();b++) {
										ComplexType sum = velocity(gf,i,j,a,b)*
										       hopping;
										w(a,b) += sum;
									}
								}
							}
						}
					}
				}
			}
		}

		template<typename GreenFunctionType>
		ComplexType velocity(const GreenFunctionType& gf,
				SizeType i,
				SizeType j,
				SizeType a,
				SizeType b) const
		{
			return std::conj(gf.matrix(i,a)) * gf.matrix(j,b) -
					std::conj(gf.matrix(j,a)) * gf.matrix(i,b);
		}

		const EngineParamsType& engineParams_;
		ParametersModelType mp_;
		const GeometryType& geometry_;
		DynVarsType dynVars_;
		SizeType hilbertSize_;
		ProgressIndicatorType progress_;
		SpinOperationsType spinOperations_;
		ThreeOrbitalTermsType threeOrbitalTerms_;
		PnictidesMultiOrbitalsObsStoredType PnictidesMultiOrbitalsObsStored_;
	}; // PnictidesMultiOrbitals

	template<typename EngineParamsType,
	         typename GeometryType>
	std::ostream& operator<<(std::ostream& os,
	                         const PnictidesMultiOrbitals<
	                           EngineParamsType,
	                           GeometryType>& model)
	{
		os<<"ModelParameters\n";
		os<<model.mp_;
		return os;
	}
} // namespace Spf

/*@}*/
#endif // PNICTIDES_MULTI_ORBS_H
