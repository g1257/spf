
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

namespace Spf {
	template<typename EngineParametersType,typename ModelType,typename RngType>
	class AlgorithmDiag {
		static const bool DO_GLAUBER = true;

	public:	
		typedef typename EngineParametersType::FieldType RealType;
		typedef std::complex<RealType> ComplexType;
		typedef PsimagLite::Matrix<ComplexType> MatrixType;

		AlgorithmDiag(const EngineParametersType& engineParams,ModelType& model)
			: engineParams_(engineParams),model_(model),
					eigNew_(model.hilbertSize()),eigOld_(model.hilbertSize()),
					hilbertSize_(model_.hilbertSize()),
					matrixNew_(hilbertSize_,hilbertSize_),
					matrixOld_(hilbertSize_,hilbertSize_)
		{
		}

		void init()
		{
			model_.createHamiltonian(matrixOld_,ModelType::OLDFIELDS);
			diag(matrixOld_,eigOld_,'N');
			sort(eigOld_.begin(), eigOld_.end(), std::less<RealType>());
		}

		ModelType& model() { return model_; } // should be const

		size_t hilbertSize() const { return hilbertSize_; }

		bool isAccepted(size_t i,RngType& rng)
		{
			RealType dsDirect = model_.deltaDirect(i);
				
			//RealType oldmu=engineParams_.mu;
			
			model_.createHamiltonian(matrixNew_,ModelType::NEWFIELDS);
			diagonalize(matrixNew_,eigNew_,'N');
			//testMatrix();
			//testEigs();
			
			if (engineParams_.carriers>0) model_.adjustChemPot(eigNew_); //changes engineParams_.mu
			RealType integrationMeasure = model_.integrationMeasure(i);
				
			RealType dS = computeDeltaAction(dsDirect,integrationMeasure);
			return metropolisOrGlauber(dS,rng);
		}

		void accept(size_t i)
		{
			model_.accept(i);
			eigOld_ = eigNew_;
		}

		void prepare()
		{
			diagonalize(matrixNew_,eigNew_,'V',ModelType::OLDFIELDS);
		}
		
		const ComplexType& matrix(size_t lambda1,size_t lambda2) const
		{
			return matrixNew_(lambda1,lambda2);
		}
		
		const RealType& e(size_t i) const
		{
			return eigNew_[i];
		}

		void diagonalize(
				MatrixType& matrix,
				std::vector<RealType>& eigs,
				char jobz='N',
				size_t fields=ModelType::NEWFIELDS) const
		{
			model_.createHamiltonian(matrix,fields);
			diag(matrix,eigs,jobz);
			if (jobz!='V')
				std::sort(eigs.begin(), eigs.end(), std::less<RealType>());
		}
		
		void printMatrix(size_t mode) const
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
		friend std::ostream& operator<<
			(std::ostream& os,const AlgorithmDiag<EngineParametersType2,
					ModelType2,RandomNumberGeneratorType2>& a);

	private:
		bool computeDeltaAction(RealType dsDirect,
		                  RealType integrationMeasure) const
		{
			RealType mu=engineParams_.mu;
			RealType beta = engineParams_.beta;
			RealType X =1.0;
			
			for (size_t i=0;i<eigNew_.size();i++) {
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
			return log(X)-beta*dsDirect;
		}
		
		bool metropolisOrGlauber(const RealType& dS2,RngType& rng) const
		{
			RealType dS = dS2;
			if (DO_GLAUBER) {
				if (dS<0) {
					dS=exp(dS)/(1.0+exp(dS));
				} else {
					dS=1.0/(1.0+exp(-dS));
				}
				return (dS>rng());
			}
			// METROPOLIS PROPER 
			return (dS > 0.0 || rng() < exp (dS));
		}
		
		void testEigs() const
		{
			RealType eps = 1e-6;
			for (size_t i=0;i<eigOld_.size();i++) {
				if (fabs(eigOld_[i]-eigNew_[i])>eps) return;
			}
			throw std::runtime_error("Eigs are equal!!\n");
		}
		
		void testMatrix() const
		{
			RealType eps = 1e-6;
			for (size_t i=0;i<matrixOld_.n_row();i++) {
				for (size_t j=0;j<matrixOld_.n_col();j++) {
					if (fabs(real(matrixOld_(i,j)-matrixNew_(i,j)))>eps &&
					fabs(imag(matrixOld_(i,j)-matrixNew_(i,j)))>eps) return;
				}
			}
			throw std::runtime_error("Matrix are equal!!\n");
		}

		const EngineParametersType& engineParams_;
		ModelType& model_;
		std::vector<RealType> eigNew_,eigOld_;
		size_t hilbertSize_;
		MatrixType matrixNew_,matrixOld_;
	}; // AlgorithmDiag
	
	template<typename EngineParametersType,typename ModelType,
		typename RandomNumberGeneratorType>
	std::ostream& operator<<(std::ostream& os,const AlgorithmDiag<
			EngineParametersType,ModelType,RandomNumberGeneratorType>& a)
	{
		
		typedef typename EngineParametersType::FieldType RealType;
		std::vector<RealType> eigNew(a.hilbertSize_);
		PsimagLite::Matrix<std::complex<RealType> > matrix(a.hilbertSize_,a.hilbertSize_);
		a.diagonalize(matrix,eigNew,'V',ModelType::OLDFIELDS);
		os<<"Eigenvalues\n";
		os<<eigNew;
		return os;
	}
} // namespace Spf

/*@}*/
#endif // ALGORITHM_DIAG_H
