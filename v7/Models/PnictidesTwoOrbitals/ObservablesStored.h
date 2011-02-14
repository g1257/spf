
/** \ingroup SPF */
/*@{*/

/*! \file ObservablesStored.h
 *
 *  ObservablesStored for pnictidesTwoOrbitals model
 *
 */

#ifndef OBSERVABLES_STORED_H
#define OBSERVABLES_STORED_H
#include "Vector.h"

namespace Spf {
	template<typename SpinOperationsType,typename ComplexType>
	class ObservablesStored {
		
		typedef typename SpinOperationsType::DynVarsType DynVarsType;
		typedef typename DynVarsType::FieldType FieldType;
		typedef PsimagLite::Vector<FieldType> VectorType;
		typedef PsimagLite::Vector<ComplexType> ComplexVectorType;
		typedef typename SpinOperationsType::GeometryType GeometryType;
		enum {DIRECTION_X,DIRECTION_Y,DIRECTION_Z};
		enum {SPIN_UP,SPIN_DOWN};
		static size_t const DIRECTIONS  = 3;

	public:
		static size_t const ORBITALS = 2;
		
		ObservablesStored(
				SpinOperationsType& spinOperations,
				const GeometryType& geometry,
				size_t dof) :
			spinOperations_(spinOperations),
			geometry_(geometry),
			dof_(dof),
			cc_(geometry.volume(),0),
			lc_(dof*geometry.volume(),0),
			chargeCor_(geometry.volume(),0),
			mc_(geometry.volume(),0),
			tc_(geometry.volume(),0),
			counter_(0)
		{}
				
		template<typename GreenFunctionType>
		void operator()(const DynVarsType& spins,
				GreenFunctionType& greenFunction)
		{
			spinOperations_.classicalCorrelations(cc_,spins);
			greenFunction.localCharge(lc_);
			chargeCorrelation(chargeCor_,greenFunction);
			mCorrelation(mc_,spins,greenFunction);
			//tCorrelation(tc_,greenFunction);
			counter_++;
		}
		
		void finalize(std::ostream& fout)
		{
			divideAndPrint(fout,cc_,"#ClassicalCorrelations:");
			divideAndPrint(fout,lc_,"#LocalCharge:");
			divideAndPrint(fout,chargeCor_,"#ChargeCorrelations:");
			divideAndPrint(fout,mc_,"#MCorrelations:");
		}

	private:

		template<typename GreenFunctionType>
		void mCorrelation(
				PsimagLite::Vector<ComplexType>& cc,
				const DynVarsType& spins,
				GreenFunctionType& greenFunction)
		{
			psimag::Matrix<ComplexType> v(spins.size,DIRECTIONS);
			calcMagElectrons(v,greenFunction);
			psimag::Matrix<FieldType> mi(spins.size,DIRECTIONS);
			calcMagSpins(mi,spins);
			for (size_t i=0;i<v.n_row();i++)
				for (size_t dir=0;dir<DIRECTIONS;dir++)
					v(i,dir) += mi(i,dir);
			correlation(mc_,v,greenFunction);
		}

		template<typename GreenFunctionType>
		void tCorrelation(
				PsimagLite::Vector<FieldType>& cc,
				GreenFunctionType& greenFunction)
		{

			psimag::Matrix<FieldType> ti(cc.size(),DIRECTIONS);
			calcTi(ti,greenFunction);
			correlation(tc_,ti,greenFunction);
		}

		// C_x = \sum_i <M_i M_{i+x}>, where M_i is passed in
		template<typename GreenFunctionType>
		void correlation(
				PsimagLite::Vector<ComplexType>& cc,
				const psimag::Matrix<ComplexType>& m,
				GreenFunctionType& greenFunction)
		{
			for (size_t x=0;x<cc.size();x++) {
				for (size_t i=0;i<cc.size();i++) {
					size_t j = geometry_.add(i,x);
					for (size_t dir=0;dir<DIRECTIONS;dir++)
						cc[x] += m(i,dir) * m(j,dir);
				}
			}
		}

		// Magnetization vector of the electrons
		template<typename GreenFunctionType>
		void calcMagElectrons(
				psimag::Matrix<ComplexType>& me,
				GreenFunctionType& greenFunction)
		{

			ComplexType sqrtOfMinus1 = ComplexType(0,1);
			size_t volume = me.n_row();
			for (size_t i=0;i<volume;i++) {
				for (size_t orb=0;orb<ORBITALS;orb++) {
					size_t x = i+(orb+SPIN_UP*ORBITALS)*volume;
					size_t y = i+(orb+SPIN_DOWN*ORBITALS)*volume;

					me(i,DIRECTION_X) -= 0.5*(
							greenFunction(x,y) + greenFunction(y,x));
					me(i,DIRECTION_Y) -= 0.5*sqrtOfMinus1*(
							greenFunction(x,y) - greenFunction(y,x));
					me(i,DIRECTION_Z) -= 0.5*
							(greenFunction(x,x) - greenFunction(y,y));
				}
			}
		}

		// Magnetization of the classical spins
		void calcMagSpins(
				psimag::Matrix<FieldType>& ms,
				const DynVarsType& spins)
		{
			size_t volume = ms.n_row();
			for (size_t i=0;i<volume;i++) {
				FieldType theta = spins.theta[i];
				FieldType phi = spins.phi[i];
				ms(i,DIRECTION_X) = sin(theta) * cos(phi);
				ms(i,DIRECTION_Y) = sin(theta) * sin(phi);
				ms(i,DIRECTION_X) = cos(theta);
			}
		}
		// C_x = \sum_i <n_i n_{i+x}>, where n_i = \sum_{dof} n_{i dof}  
		template<typename GreenFunctionType>
		void chargeCorrelation(PsimagLite::Vector<FieldType>& cc,
				       GreenFunctionType& greenFunction)
		{
			size_t volume = geometry_.volume();
			//allChargeCorrelation(greenFunction); // for debugging only, comment out for production
			for (size_t gamma=0;gamma<dof_;gamma++) {
				for (size_t gamma2=0;gamma2<dof_;gamma2++) {
					PsimagLite::Vector<FieldType> tmpV(dof_*dof_*volume,0);
					chargeCorrelation(tmpV,gamma,gamma2,greenFunction);
					for (size_t x=0;x<volume;x++) {
						cc[x] += tmpV[x + gamma*volume+gamma2*volume*dof_];
					}
				}
			}
		}
		
		template<typename GreenFunctionType>
		void chargeCorrelation(PsimagLite::Vector<FieldType>& cc,
				       size_t gamma,size_t gamma2,GreenFunctionType& greenFunction)
		{
			size_t volume = geometry_.volume();
			size_t dof = dof_; 
        		// Do for all x's
			FieldType slip = 0;
        		for (size_t x=0;x<volume;x++) {     // Sum over all i's (or w's)
                		ComplexType tmp=0;
				//ComplexType test= ComplexType(0.,0.);
                		for (size_t w=0;w<volume;w++) {
					size_t v=geometry_.add(x,w);
                                        tmp += chargeCorrelation(w,gamma,v,gamma2,greenFunction);
					
					//test += (1.0-greenFunction(w+gamma*volume,
					//			w+gamma*volume));
				}
				slip += imag(tmp);
                		cc[x+gamma*volume+gamma2*volume*dof] += real(tmp)/volume;
			}
			
			if (fabs(slip)>1e-6) std::cerr<<"slipping="<<slip<<"\n";
		}
		template<typename GreenFunctionType>
		void allChargeCorrelation(GreenFunctionType& greenFunction)
		{
			size_t volume = geometry_.volume();
			psimag::Matrix<FieldType> m(volume*2,volume*2);
			for (size_t i=0;i<volume;i++) {
				for (size_t gamma=0;gamma<4;gamma++) {
					size_t orb1 = (gamma & 1);
					for (size_t j=0;j<volume;j++) {
						for (size_t gamma2=0;gamma2<4;gamma2++) {
							size_t orb2 = (gamma2 & 1);
							m(i+orb1*volume,j+orb2*volume) += 
								real(chargeCorrelation(i,gamma,j,gamma2,greenFunction));
						}
					}
				}
			}
			std::cerr<<m;
			throw std::runtime_error("testing\n");
							
		}
		
		template<typename GreenFunctionType>
		ComplexType chargeCorrelation(size_t w, size_t gamma,size_t v, size_t gamma2,
			GreenFunctionType& greenFunction)
		{
			size_t volume = geometry_.volume();
			ComplexType tmp = 0;
			if (v==w && gamma==gamma2) {
				tmp += (1.0-greenFunction(w+gamma*volume,w+gamma*volume));
			} else {
				tmp +=(1.0-greenFunction(w+gamma*volume,w+gamma*volume))*
					(1.0-greenFunction(v+gamma2*volume,v+gamma2*volume));
				tmp -= greenFunction(v+gamma2*volume,w+gamma*volume)*
					greenFunction(w+gamma*volume,v+gamma2*volume);
			}
			return tmp;
		}
		
		template<typename SomeVectorType>
		void divideAndPrint(
				std::ostream& fout,
				SomeVectorType& v,
				const std::string& label)
		{
			v /= counter_;
			fout<<label<<"\n";
			fout<<v;
		}

		SpinOperationsType& spinOperations_;
		const GeometryType& geometry_;
		size_t dof_;
		VectorType cc_;
		VectorType lc_;
		VectorType chargeCor_;
		ComplexVectorType mc_;
		ComplexVectorType tc_;
		size_t counter_;

	}; // ObservablesStored
	
} // namespace Spf

/*@}*/
#endif // OBSERVABLES_STORED_H
