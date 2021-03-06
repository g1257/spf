
/** \ingroup SPF */
/*@{*/

/*! \file DmsMultiOrbitalObsStored.h
 *
 *  DmsMultiOrbitalObsStored for DmsMultiOrbital model
 *
 */

#ifndef OBS_STORED_DMS_MULTI_ORB_H
#define OBS_STORED_DMS_MULTI_ORB_H
#include "Histogram.h"
#include "Vector.h" // in PsimagLite
#include "Matrix.h" // in PsimagLite
#include "Fermi.h"  // in PsimagLite

namespace Spf {
template<typename SpinOperationsType,
         typename ComplexType,
         typename ParametersModelType,
         typename EngineParamsType>
class DmsMultiOrbitalObsStored {

	typedef typename SpinOperationsType::DynVarsType DynVarsType;
	typedef typename DynVarsType::FieldType RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename SpinOperationsType::GeometryType GeometryType;
	typedef PsimagLite::Matrix<RealType> MatrixType;
	typedef Histogram<RealType,ComplexType> HistogramComplexType;
	typedef Histogram<RealType,RealType> HistogramRealType;

	enum {DIRECTION_X,DIRECTION_Y,DIRECTION_Z};
	static SizeType const DIRECTIONS  = 3;

public:

	DmsMultiOrbitalObsStored(
	        SpinOperationsType& spinOperations,
	        const GeometryType& geometry,
	        const ParametersModelType& mp,
	        const EngineParamsType& pe)
	    : spinOperations_(spinOperations),
	      geometry_(geometry),
	      mp_(mp),
	      pe_(pe),
	      cs_(geometry.volume(),DIRECTIONS),
	      arw_(0),
	      optical_(0),
	      counter_(0)
	{
		if (pe_.options.find("akw")!=PsimagLite::String::npos)
			arw_.resize(geometry.volume()*mp_.orbitals*2,
			            new HistogramComplexType(mp.histogramParams[0],
			            mp.histogramParams[1],
			        SizeType(mp.histogramParams[2])));

		if (pe_.options.find("optical")!=PsimagLite::String::npos)
			optical_ = new HistogramRealType(mp.histogramParams[0],
			        mp.histogramParams[1],
			        SizeType(mp.histogramParams[2]));
	}

	~DmsMultiOrbitalObsStored()
	{
		for (SizeType i = 0; i < arw_.size(); ++i)
			if (arw_[i]) delete arw_[i];

		if (optical_) delete optical_;
	}

	template<typename GreenFunctionType>
	void operator()(const DynVarsType& spins,
	                GreenFunctionType& greenFunction)
	{
		PsimagLite::Matrix<FieldType> mi(spins.size,DIRECTIONS);
		calcMagSpins(mi,spins);
		correlation(cs_,mi);

		if (pe_.options.find("akw")!=PsimagLite::String::npos)
			accAkw(greenFunction);
		if (pe_.options.find("optical")!=PsimagLite::String::npos)
			accOptical(greenFunction);
		counter_++;
	}

	template<typename SomeOutputType>
	void finalize(SomeOutputType& fout)
	{
		divideAndPrint(fout,cs_,"#ClassicalSpinCorrelations");

		if (pe_.options.find("akw")!=PsimagLite::String::npos)
			reduce(arw_);
		if (pe_.options.find("optical")!=PsimagLite::String::npos)
			optical_->reduce();

		if (!PsimagLite::Concurrency::root()) return;

		if (pe_.options.find("akw")!=PsimagLite::String::npos)
			divideAndPrint(fout,arw_,"#Arw:");
		if (pe_.options.find("optical")!=PsimagLite::String::npos)
			divideAndPrint(fout,optical_,"#Optical:");
	}

private:

	void calcMagSpins(PsimagLite::Matrix<FieldType>& ms,
	                  const DynVarsType& spins)
	{
		SizeType volume = ms.n_row();
		for (SizeType i=0;i<volume;i++) {
			FieldType theta = spins.theta[i];
			FieldType phi = spins.phi[i];
			FieldType m = (mp_.modulus[i]) ? 1.0 : 0.0;
			ms(i,DIRECTION_X) = m*sin(theta) * cos(phi);
			ms(i,DIRECTION_Y) = m*sin(theta) * sin(phi);
			ms(i,DIRECTION_Z) = m*cos(theta);
		}
	}

	void correlation(MatrixType& cc,
	                 const PsimagLite::Matrix<FieldType>& m)
	{
		for (SizeType x=0;x<cc.n_row();x++) {
			for (SizeType i=0;i<cc.n_row();i++) {
				SizeType j = geometry_.add(i,x);
				for (SizeType dir=0;dir<DIRECTIONS;dir++)
					cc(x,dir) += std::real(m(i,dir) * m(j,dir));
			}
		}
	}

	void reduce(typename PsimagLite::Vector<HistogramComplexType*>::Type& h)
	{
		for (SizeType i=0;i<h.size();i++) h[i]->reduce();
	}

	//! A(r+gamma*N,omega) will contain A(r,omega)_\gamma
	template<typename GreenFunctionType>
	void accAkw(const GreenFunctionType& gf)
	{
		SizeType n = geometry_.volume();
		SizeType dof = 2*mp_.orbitals;

		for (SizeType r=0;r<n;r++) {
			for (SizeType l=0;l<gf.hilbertSize();l++) {
				for (SizeType gamma=0;gamma<dof;gamma++) {
					ComplexType temp = 0.0;
					for (SizeType i=0;i<n;i++) {
						SizeType j=geometry_.add(i,r);
						temp += conj(gf.matrix(i+gamma*n,l))*
						        gf.matrix(j+gamma*n,l);
						arw_[r+gamma*n]->add(gf.e(l),temp);
					}
				}
			}
		}
	}

	template<typename GreenFunctionType>
	void accOptical(const GreenFunctionType& gf)
	{
		SizeType n = geometry_.volume();
		SizeType dof  = 2*mp_.orbitals;
		//int aindex = 0*(ether.linSize/4); //(0,0,0)
		//int aindex = 1*(ether.linSize/4); //(0,a/2,a/2)
		//int aindex = 2*(ether.linSize/4); //(a/2,0,a/2)
		SizeType aindex = 3*(n/4); //(a/2,a/2,0)

		// some checking
		if (geometry_.name()!="fcc" || n % 4!=0) {
			PsimagLite::String s = "accOptical: "+ PsimagLite::String(__FILE__) +
			        " " + ttos(__LINE__) + "\n";
			throw PsimagLite::RuntimeError(s.c_str());
		}

		SizeType hilbertSize = gf.hilbertSize();
		RealType beta = pe_.beta;

		for (SizeType l=0;l<hilbertSize;l++) {
			RealType e1 = gf.e(l)-pe_.mu;
			for (SizeType l2=0;l2<hilbertSize;l2++) {
				RealType e2 = gf.e(l2)-pe_.mu;
				if (e1-e2<1e-3) continue;
				ComplexType temp = 0.0;
				for (SizeType i=0;i<n;i++) {
					SizeType j = geometry_.add(i,aindex); // j = i+avector
					SizeType dir = geometry_.scalarDirection(i,j);
					for (SizeType gamma=0;gamma<dof;gamma++) {
						for (SizeType gamma2=0;gamma2<dof;gamma2++) {
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
				optical_->add(e1-e2,temp2);
			}
		}
	}

	template<typename SomeOutputType>
	void divideAndPrint(SomeOutputType& fout,
	                    VectorRealType& v,
	                    const PsimagLite::String& label)
	{
		v /= counter_;
		if (!PsimagLite::Concurrency::root()) return;
		fout<<label<<"\n";
		fout<<v;
	}

	template<typename SomeOutputType>
	void divideAndPrint(SomeOutputType& fout,
	                    MatrixType& m,
	                    const PsimagLite::String& label)
	{
		VectorRealType v(m.n_row(),0.0);
		for (SizeType dir=0;dir<m.n_col();dir++) {
			for (SizeType i=0;i<m.n_row();i++) v[i] =  m(i,dir);
			PsimagLite::String newlabel = label+ttos(dir);
			divideAndPrint(fout,v,newlabel);
		}
	}

	template<typename SomeOutputType>
	void divideAndPrint(SomeOutputType& fout,
	                    typename PsimagLite::Vector<HistogramComplexType*>::Type& h,
	                    const PsimagLite::String& label)
	{
		if (h.size()==0) return;

		for (SizeType i=0;i<h.size();i++) {
			h[i]->divide(counter_*h[i]->xWidth());
		}
		fout<<label<<"\n";
		for (SizeType i=0;i<h[0]->size();i++)  {
			fout<<h[0]->x(i)<<" ";
			for (SizeType j=0;j<h.size();j++)
				fout<<std::real(h[j]->y(i))<<" ";
			fout<<"\n";
		}
		fout<<"% Imaginary part below\n";
		for (SizeType i=0;i<h[0]->size();i++)  {
			fout<<h[0]->x(i)<<" ";
			for (SizeType j=0;j<h.size();j++)
				fout<<std::imag(h[j]->y(i))<<" ";
			fout<<"\n";
		}
	}

	template<typename SomeOutputType>
	void divideAndPrint(SomeOutputType& fout,
	                    HistogramRealType* h,
	                    const PsimagLite::String& label)
	{

		h->divide(counter_*h->xWidth());
		fout<<label<<"\n";
		for (SizeType i=0;i<h->size();i++)
			fout<<h->x(i)<<" "<<h->y(i)<<"\n";
	}

	SpinOperationsType& spinOperations_;
	const GeometryType& geometry_;
	const ParametersModelType& mp_;
	const EngineParamsType& pe_;
	MatrixType cs_;
	typename PsimagLite::Vector<HistogramComplexType*>::Type arw_;
	HistogramRealType* optical_;
	SizeType counter_;
}; // DmsMultiOrbitalObsStored

} // namespace Spf

/*@}*/
#endif // OBS_STORED_DMS_MULTI_ORB_H

