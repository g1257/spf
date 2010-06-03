
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
		typedef typename SpinOperationsType::GeometryType GeometryType;
		
	public:
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
			counter_(0)
		{}
				
		template<typename GreenFunctionType>
		void operator()(const DynVarsType& spins,GreenFunctionType& greenFunction)
		{
			spinOperations_.classicalCorrelations(cc_,spins);
			greenFunction.localCharge(lc_);
			chargeCorrelation(chargeCor_,greenFunction);
			counter_++;
		}
		
		// C_x = \sum_i <n_i n_{i+x}>, where n_i = \sum_{dof} n_{i dof}  
		template<typename GreenFunctionType>
		void chargeCorrelation(PsimagLite::Vector<FieldType>& cc,
				       GreenFunctionType& greenFunction)
		{
			size_t volume = geometry_.volume();
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
			size_t dof = dof_; // the 2 comes because of the spin
			
        		// Do for all x's
        		for (size_t x=0;x<volume;x++) {     // Sum over all i's (or w's)
                		ComplexType tmp=0;
                		for (size_t w=0;w<volume;w++) {
					size_t v=geometry_.add(x,w);
                                        if (x==0 && gamma==gamma2) {
                                                tmp += (1.0-greenFunction(w+gamma*volume,
								w+gamma*volume));
                                        } else {
						tmp +=(1.0-greenFunction(w+gamma*volume,
								w+gamma*volume))*
						(1.0-greenFunction(v+gamma2*volume,
									v+gamma2*volume));
                                                tmp -= greenFunction(v+gamma2*volume,
								w+gamma*volume)*
						greenFunction(w+gamma*volume,v+gamma2*volume);
					}
				}
                		cc[x+gamma*volume+gamma2*volume*dof] += real(tmp)/volume;
			}
		}
		
		void finalize(std::ostream& fout)
		{
			divideAndPrint(fout,cc_,"#ClassicalCorrelations:");
			divideAndPrint(fout,lc_,"#LocalCharge:");
			divideAndPrint(fout,chargeCor_,"#ChargeCorrelations:");
		}
		
	private:
		SpinOperationsType& spinOperations_;
		const GeometryType& geometry_;
		size_t dof_;
		VectorType cc_;
		VectorType lc_;
		VectorType chargeCor_;
		size_t counter_;
		
		void divideAndPrint(std::ostream& fout,VectorType& v,const std::string& label)
		{
			v /= counter_;
			fout<<label<<"\n";
			fout<<v;
		}

	}; // ObservablesStored
	
} // namespace Spf

/*@}*/
#endif // OBSERVABLES_STORED_H
