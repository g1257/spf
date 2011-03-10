
/** \ingroup SPF */
/*@{*/

/*! \file ObservablesStored.h
 *
 *  ObservablesStored for DmsMultiOrbital model
 *
 */

#ifndef OBS_STORED_DMS_MULTI_ORB_H
#define OBS_STORED_DMS_MULTI_ORB_H
#include "Vector.h"

namespace Spf {
	template<typename SpinOperationsType,typename ComplexType>
	class ObservablesStored {
		
		typedef typename SpinOperationsType::DynVarsType DynVarsType;
		typedef typename DynVarsType::FieldType FieldType;
		typedef PsimagLite::Vector<FieldType> VectorType;
		typedef PsimagLite::Vector<ComplexType> ComplexVectorType;
		typedef typename SpinOperationsType::GeometryType GeometryType;
		typedef PsimagLite::Matrix<FieldType> MatrixType;

		//enum {DIRECTION_X,DIRECTION_Y,DIRECTION_Z};
		//enum {ORBITAL_XZ,ORBITAL_YZ};
		//enum {SPIN_UP,SPIN_DOWN};
		//static size_t const DIRECTIONS  = 3;

	public:
		static size_t const ORBITALS = 3;
		
		ObservablesStored(
				SpinOperationsType& spinOperations,
				const GeometryType& geometry,
				size_t dof) :
			spinOperations_(spinOperations),
			geometry_(geometry),
			dof_(dof),
//			lc_(dof*geometry.volume(),0),
//			chargeCor_(geometry.volume(),0),
//			mc_(geometry.volume(),DIRECTIONS),
//			tc_(geometry.volume(),DIRECTIONS),
//			cs_(geometry.volume(),DIRECTIONS),
//			qs_(geometry.volume(),DIRECTIONS),
			counter_(0)
		{}
				
		template<typename GreenFunctionType>
		void operator()(const DynVarsType& spins,
				GreenFunctionType& greenFunction)
		{
//			greenFunction.localCharge(lc_);
//			chargeCorrelation(chargeCor_,greenFunction);
//			psimag::Matrix<FieldType> mi(spins.size,DIRECTIONS);
//			calcMagSpins(mi,spins);
//			psimag::Matrix<ComplexType> qi(spins.size,DIRECTIONS);
//			calcMagElectrons(qi,greenFunction);
//			mCorrelation(mi,qi,greenFunction);
//			tCorrelation(greenFunction);
//			correlation(cs_,mi,greenFunction);
//			correlation(qs_,qi,greenFunction);
			counter_++;
		}
		
		void finalize(std::ostream& fout)
		{
//			divideAndPrint(fout,lc_,"#LocalCharge:");
//			divideAndPrint(fout,chargeCor_,"#ChargeCorrelations:");
//			divideAndPrint(fout,mc_,"#MCorrelations");
//			divideAndPrint(fout,tc_,"#TCorrelations");
//			divideAndPrint(fout,cs_,"#ClassicalSpinCorrelations");
//			divideAndPrint(fout,qs_,"#ItinerantSpinCorrelations");
		}

	private:

		
		void divideAndPrint(
				std::ostream& fout,
				VectorType& v,
				const std::string& label)
		{
			v /= counter_;
			fout<<label<<"\n";
			fout<<v;
		}

		void divideAndPrint(
				std::ostream& fout,
				MatrixType& m,
				const std::string& label)
		{
			VectorType v(m.n_row(),0);
			for (size_t dir=0;dir<m.n_col();dir++) {
				for (size_t i=0;i<m.n_row();i++) v[i] =  m(i,dir);
				std::string newlabel = label+utils::ttos(dir);
				divideAndPrint(fout,v,newlabel);
			}
		}

		SpinOperationsType& spinOperations_;
		const GeometryType& geometry_;
		size_t dof_;
//		VectorType lc_;
//		VectorType chargeCor_;
//		MatrixType mc_;
//		MatrixType tc_;
//		MatrixType cs_;
//		MatrixType qs_;
		size_t counter_;

	}; // ObservablesStored
	
} // namespace Spf

/*@}*/
#endif // OBS_STORED_DMS_MULTI_ORB_H
