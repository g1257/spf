
/** \ingroup SPF */
/*@{*/

/*! \file ObservablesStored.h
 *
 *  ObservablesStored for DmsMultiOrbital model
 *
 */

#ifndef OBS_STORED_DMS_MULTI_ORB_H
#define OBS_STORED_DMS_MULTI_ORB_H
#include "Histogram.h"
#include "Vector.h" // in PsimagLite
#include "Matrix.h" // in PsimagLite
#include "Fermi.h"  // in PsimagLite

namespace Spf {
	template<typename SpinOperationsType,typename ComplexType,
	         typename ParametersModelType,typename EngineParamsType,
	         typename ConcurrencyType>
	class ObservablesStored {
		
		typedef typename SpinOperationsType::DynVarsType DynVarsType;
		typedef typename DynVarsType::FieldType RealType;
//		typedef PsimagLite::Vector<FieldType> VectorType;
		typedef PsimagLite::Vector<ComplexType> ComplexVectorType;
		typedef typename SpinOperationsType::GeometryType GeometryType;
		typedef PsimagLite::Matrix<RealType> MatrixType;

		typedef Histogram<RealType,ComplexType> HistogramComplexType;
		typedef Histogram<RealType,RealType> HistogramRealType;
		typedef typename ConcurrencyType::CommType CommType;

	public:
		static size_t const ORBITALS = 3;
		
		ObservablesStored(
				SpinOperationsType& spinOperations,
				const GeometryType& geometry,
				const ParametersModelType& mp,
				const EngineParamsType& pe,
				ConcurrencyType& concurrency)
 		: spinOperations_(spinOperations),
		 geometry_(geometry),
		 mp_(mp),
		 pe_(pe),
		 concurrency_(concurrency),
		 arw_(geometry.volume()*ORBITALS*2,HistogramComplexType(
		       mp.histogramParams[0],mp.histogramParams[1],
		       size_t(mp.histogramParams[2]))),
		       optical_(mp.histogramParams[0],mp.histogramParams[1],
		       size_t(mp.histogramParams[2])),
		 counter_(0)
		{}

		template<typename GreenFunctionType>
		void operator()(const DynVarsType& spins,
				GreenFunctionType& greenFunction)
		{
			if (pe_.options.find("akw")!=std::string::npos)
				accAkw(greenFunction);
			if (pe_.options.find("optical")!=std::string::npos)
				accOptical(greenFunction);
			counter_++;
		}
		
		template<typename SomeOutputType>
		void finalize(SomeOutputType& fout,CommType comm)
		{
			reduce(arw_,comm);
			optical_.reduce(concurrency_,comm);
			if (!concurrency_.root(comm)) return;
			divideAndPrint(fout,arw_,"#Arw:");
			divideAndPrint(fout,optical_,"#Optical:");
		}

	private:

		void reduce(PsimagLite::Vector<HistogramComplexType>::Type& h,CommType comm)
		{
			for (size_t i=0;i<h.size();i++) h[i].reduce(concurrency_,comm);
		}

		//! A(r+gamma*N,omega) will contain A(r,omega)_\gamma
		template<typename GreenFunctionType>
		void accAkw(const GreenFunctionType& gf)
		{
			size_t n = geometry_.volume();
			size_t dof = 2*ORBITALS;

			for (size_t r=0;r<n;r++) {
				for (size_t l=0;l<gf.hilbertSize();l++) {
					for (size_t gamma=0;gamma<dof;gamma++) {
						ComplexType temp = 0.0;
						for (size_t i=0;i<n;i++) {
							size_t j=geometry_.add(i,r);
							temp += conj(gf.matrix(i+gamma*n,l))*
									gf.matrix(j+gamma*n,l);
							arw_[r+gamma*n].add(gf.e(l),temp);
						}
					}
				}
			}
		}
		
		template<typename GreenFunctionType>
		void accOptical(const GreenFunctionType& gf)
		{
			size_t n = geometry_.volume();
			size_t dof  = 2*ORBITALS;
			//int aindex = 0*(ether.linSize/4); //(0,0,0)
			//int aindex = 1*(ether.linSize/4); //(0,a/2,a/2)
			//int aindex = 2*(ether.linSize/4); //(a/2,0,a/2)
			size_t aindex = 3*(n/4); //(a/2,a/2,0)


			// some checking
			if (geometry_.name()!="fcc" || n % 4!=0) {
				std::string s = "accOptical: "+ std::string(__FILE__) +
						" " + ttos(__LINE__) + "\n";
				throw std::runtime_error(s.c_str());
			}

			size_t hilbertSize = gf.hilbertSize();
			RealType beta = pe_.beta;

			for (size_t l=0;l<hilbertSize;l++) {
				RealType e1 = gf.e(l)-pe_.mu;
				for (size_t l2=0;l2<hilbertSize;l2++) {
					RealType e2 = gf.e(l2)-pe_.mu;
					if (e1-e2<1e-3) continue;
					ComplexType temp = 0.0;
					for (size_t i=0;i<n;i++) {
						size_t j = geometry_.add(i,aindex); // j = i+avector
						size_t dir = geometry_.scalarDirection(i,j);
						for (size_t gamma=0;gamma<dof;gamma++) {
							for (size_t gamma2=0;gamma2<dof;gamma2++) {
								ComplexType thop =
									mp_.hoppings[gamma+gamma2*dof+dir*dof*dof];
								temp += thop * conj(gf.matrix(i+gamma*n,l))*
										gf.matrix(j+gamma2*n,l2);
								temp -= thop * conj(gf.matrix(j+gamma*n,l))*
										gf.matrix(i+gamma2*n,l2);
							}
						}
					}
					RealType temp2 = real(temp)*real(temp)+imag(temp)*imag(temp);

					temp2 = temp2 *(PsimagLite::fermi(beta*e2)-
							PsimagLite::fermi(beta*e1))/(e1-e2);
					//temp = temp * fermi(beta*e2) * fermi(-beta*e1);
					//temp = temp * (1.0 - exp(-beta*(e1-e2)))/(e1-e2);
					optical_.add(e1-e2,temp2);
				}
			}
		}

		template<typename SomeOutputType>
		void divideAndPrint(SomeOutputType& fout,
		                    PsimagLite::Vector<HistogramComplexType>::Type& h,
		                    const std::string& label)
		{
			if (h.size()==0) return;

			for (size_t i=0;i<h.size();i++) {
				h[i].divide(counter_*h[i].xWidth());
			}
			fout<<label<<"\n";
			for (size_t i=0;i<h[0].size();i++)  {
				fout<<h[0].x(i)<<" ";
				for (size_t j=0;j<h.size();j++)
					fout<<std::real(h[j].y(i))<<" ";
				fout<<"\n";
			}
			fout<<"% Imaginary part below\n";
			for (size_t i=0;i<h[0].size();i++)  {
				fout<<h[0].x(i)<<" ";
				for (size_t j=0;j<h.size();j++)
					fout<<std::imag(h[j].y(i))<<" ";
				fout<<"\n";
			}
		}
		
		template<typename SomeOutputType>
		void divideAndPrint(SomeOutputType& fout,
		                    HistogramRealType& h,
		                    const std::string& label)
		{

			h.divide(counter_*h.xWidth());
			fout<<label<<"\n";
			for (size_t i=0;i<h.size();i++)
				fout<<h.x(i)<<" "<<h.y(i)<<"\n";
		}

		SpinOperationsType& spinOperations_;
		const GeometryType& geometry_;
		const ParametersModelType& mp_;
		const EngineParamsType& pe_;
		ConcurrencyType& concurrency_;
		PsimagLite::Vector<HistogramComplexType>::Type arw_;
		HistogramRealType optical_;
		size_t counter_;

	}; // ObservablesStored
	
} // namespace Spf

/*@}*/
#endif // OBS_STORED_DMS_MULTI_ORB_H
