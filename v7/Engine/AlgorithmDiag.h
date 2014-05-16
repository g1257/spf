
/** \ingroup SPF */
/*@{*/

/*! \file AlgorithmDiag.h
 *
 *  Diagonalization method for SPF
 *
 */
#ifndef ALGORITHM_DIAG_H
#define ALGORITHM_DIAG_H
#include <algorithm>
#include "ProgressIndicator.h" // in PsimagLite
#include "Matrix.h" // in PsimagLite
#include "Fermi.h" // in PsimagLite
#include "Complex.h" // in PsimagLite
#include "MetropolisOrGlauber.h"
#include "Adjustments.h"
#include "Sort.h"

namespace Spf {
	template<typename EngineParametersType,typename ModelType,typename RngType>
	class AlgorithmDiag {
		typedef Adjustments<EngineParametersType> AdjustmentsType;
	public:
		typedef typename EngineParametersType::RealType RealType;
		typedef typename EngineParametersType::IoInType IoInType;
		typedef typename ModelType::MatrixType MatrixType;
		typedef typename MatrixType::value_type ComplexOrRealType;
		typedef MetropolisOrGlauber<RealType,RngType> MetropolisOrGlauberType;

		AlgorithmDiag(const EngineParametersType& engineParams,
		              ModelType& model,
		              IoInType& io)
		: engineParams_(engineParams),
		  model_(model),
		  metropolisOrGlauber_(engineParams.detailedBalance),
		  adjustments_(engineParams),
		  eigNew_(model.hilbertSize()),
		  eigOld_(model.hilbertSize()),
		  hilbertSize_(model_.hilbertSize()),
		  matrixNew_(hilbertSize_,hilbertSize_),
		  matrixOld_(hilbertSize_,hilbertSize_)
		{
		}

		void init()
		{
			model_.createHamiltonian(matrixOld_,ModelType::OLDFIELDS);

			if (engineParams_.options.find("matrixBlocked")!=PsimagLite::String::npos) {
				diagBlocked(matrixOld_,eigOld_,'N');
			} else {
				diag(matrixOld_,eigOld_,'N');
			}

			sort(eigOld_.begin(), eigOld_.end(), std::less<RealType>());
		}

		ModelType& model() { return model_; } // should be const

		SizeType hilbertSize() const { return hilbertSize_; }

		template<typename OperationsType>
		bool isAccepted(SizeType i,RngType& rng,OperationsType& ops,int n)
		{
			model_.createHamiltonian(matrixNew_,ModelType::NEWFIELDS);
			diagonalize(matrixNew_,eigNew_,'N');
			//testMatrix();
			//testEigs();
			if (engineParams_.carriers>0 && needsAdjustment(i)) {
				try {
					engineParams_.mu = adjustments_.adjChemPot(eigNew_);
				} catch (std::exception& e) {
					std::cerr<<e.what()<<"\n";
				}
			}

			RealType integrationMeasure = model_.integrationMeasure(i,ops,n);
				
			RealType X = computeDeltaAction(integrationMeasure);
			X *= exp(-engineParams_.beta*model_.deltaDirect(i,ops,n));
			return metropolisOrGlauber_(X,rng);
		}

		template<typename OperationsType>
		void accept(SizeType i,OperationsType& ops)
		{
			ops.accept(i);
			eigOld_ = eigNew_;
		}

		void prepare()
		{
			diagonalize(matrixNew_,eigNew_,'V',ModelType::OLDFIELDS);
		}

		const ComplexOrRealType& matrix(SizeType lambda1,SizeType lambda2) const
		{
			return matrixNew_(lambda1,lambda2);
		}

		const RealType& e(SizeType i) const
		{
			return eigNew_[i];
		}

		void diagonalize(MatrixType& matrix,
		                 typename PsimagLite::Vector<RealType>::Type& eigs,
		                 char jobz='N',
		                 SizeType fields=ModelType::NEWFIELDS) const
		{
			model_.createHamiltonian(matrix,fields);
			if (engineParams_.options.find("matrixBlocked")!=PsimagLite::String::npos) {
				diagBlocked(matrix,eigs,jobz);
			} else {
				diag(matrix,eigs,jobz);
			}
			if (jobz!='V')
				std::sort(eigs.begin(), eigs.end(), std::less<RealType>());
		}

		void printMatrix(SizeType mode) const
		{
			if (mode==ModelType::NEWFIELDS) {
				std::cerr<<matrixNew_;
				return;
			}
			MatrixType m(matrixNew_.n_row(),matrixNew_.n_col());
			if (!isHermitian(m)) throw std::runtime_error("Problem\n");
			model_.createHamiltonian(m,ModelType::OLDFIELDS);
			std::cerr<<m;
		}

		template<typename EngineParametersType2,typename ModelType2,
			typename RandomNumberGeneratorType2>
		friend std::ostream& operator<<(std::ostream& os,
		                                const AlgorithmDiag<EngineParametersType2,
		                                ModelType2,RandomNumberGeneratorType2>& a);

	private:

		bool needsAdjustment(SizeType i) const
		{
			if (engineParams_.adjustEach==0) return true;
			SizeType x = (i % engineParams_.adjustEach);
			return (x==0);
		}

		RealType computeDeltaAction(RealType integrationMeasure) const
		{
			RealType mu=engineParams_.mu;
			RealType beta = engineParams_.beta;
			RealType X =1.0;

			for (SizeType i=0;i<eigNew_.size();i++) {
				RealType temp = 0;
				if (eigNew_[i]>mu)
					temp = (1.0+exp(-beta*(eigNew_[i]-mu)))/
						(1.0+exp(-beta*(eigOld_[i]-mu)));
				else
				temp =(1.0+exp(beta*(eigNew_[i]-mu)))/
							(exp(beta*(eigNew_[i]-mu))+
									exp(-beta*(eigOld_[i]-eigNew_[i])));
			
				X *= temp;
			}
			//std::cerr<<"Xbefore="<<X<<" ";
			//if (ether.isSet("sineupdate")) X *= integrationMeasure;
			return X;
		}

		void testEigs() const
		{
			RealType eps = 1e-6;
			for (SizeType i=0;i<eigOld_.size();i++) {
				if (fabs(eigOld_[i]-eigNew_[i])>eps) return;
			}
			throw std::runtime_error("Eigs are equal!!\n");
		}
		
		void testMatrix() const
		{
			RealType eps = 1e-6;
			for (SizeType i=0;i<matrixOld_.n_row();i++) {
				for (SizeType j=0;j<matrixOld_.n_col();j++) {
					if (fabs(real(matrixOld_(i,j)-matrixNew_(i,j)))>eps &&
					fabs(imag(matrixOld_(i,j)-matrixNew_(i,j)))>eps) return;
				}
			}
			throw std::runtime_error("Matrix are equal!!\n");
		}

		void diagBlocked(MatrixType& matrix,
		                 typename PsimagLite::Vector<RealType>::Type& eigs,
		                 char jobz) const
		{
			SizeType n = matrix.n_row();
			assert(!(n&1));
			n = SizeType(n/2);

			MatrixType m(n,n);
			for (SizeType i=0;i<n;i++)
				for (SizeType j=0;j<n;j++)
					m(i,j) = matrix(i,j);

			typename PsimagLite::Vector<RealType>::Type eigs1(n);
			diag(m,eigs1,jobz);

			MatrixType m2(n,n);
			for (SizeType i=0;i<n;i++)
				for (SizeType j=0;j<n;j++)
					m2(i,j) = matrix(i+n,j+n);

			typename PsimagLite::Vector<RealType>::Type eigs2(n);
			diag(m2,eigs2,jobz);
			combine(matrix,eigs,m,eigs1,m2,eigs2,jobz);
		}

		void combine(MatrixType& matrix,
		             typename PsimagLite::Vector<RealType>::Type& eigs,
		             MatrixType& m1,
		             typename PsimagLite::Vector<RealType>::Type& eigs1,
		             MatrixType& m2,
		             typename PsimagLite::Vector<RealType>::Type& eigs2,
		             char jobz) const
		{
			SizeType n = eigs1.size();
			assert(n==eigs2.size());
			assert(eigs.size()==2*n);
			assert(matrix.n_row()==2*n);
			assert(m1.n_row()==n);
			assert(m2.n_row()==n);

			for (SizeType i=0;i<n;i++) {
				eigs[i] = eigs1[i];
				eigs[i+n] = eigs2[i];
			}

			PsimagLite::Sort<typename PsimagLite::Vector<RealType>::Type> sort;
			PsimagLite::Vector<SizeType>::Type iperm(eigs.size());

			sort.sort(eigs,iperm);

			if (jobz!='v' && jobz!='V') return;

			MatrixType m3(2*n,2*n);
			for (SizeType i=0;i<n;i++) {
				for (SizeType j=0;j<n;j++) {
					m3(i,j) = m1(i,j);
					m3(i+n,j+n) = m2(i,j);
				}
			}
			for (SizeType i=0;i<2*n;i++) {
				for (SizeType j=0;j<2*n;j++) {
					matrix(i,j) = m3(iperm[i],iperm[j]);
				}
			}
		}

		const EngineParametersType& engineParams_;
		ModelType& model_;
		MetropolisOrGlauberType metropolisOrGlauber_;
		AdjustmentsType adjustments_;
		typename PsimagLite::Vector<RealType>::Type eigNew_,eigOld_;
		SizeType hilbertSize_;
		MatrixType matrixNew_,matrixOld_;
	}; // AlgorithmDiag
	
	template<typename EngineParametersType,typename ModelType,
		typename RandomNumberGeneratorType>
	std::ostream& operator<<(std::ostream& os,const AlgorithmDiag<
			EngineParametersType,ModelType,RandomNumberGeneratorType>& a)
	{
		
		typedef typename EngineParametersType::RealType RealType;
		typename PsimagLite::Vector<RealType>::Type eigNew(a.hilbertSize_);
		typename ModelType::MatrixType matrix(a.hilbertSize_,a.hilbertSize_);
		a.diagonalize(matrix,eigNew,'V',ModelType::OLDFIELDS);
		os<<"#Eigenvalues\n";
		os<<eigNew;
		return os;
	}
} // namespace Spf

/*@}*/
#endif // ALGORITHM_DIAG_H
