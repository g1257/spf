
/** \ingroup SPF */
/*@{*/

/*! \file AlgorithmDiag.h
 *
 *  Diagonalization method for SPF
 *
 */
#ifndef ALGORITHM_DIAG_H
#define ALGORITHM_DIAG_H
#include "Utils.h"
#include "ProgressIndicator.h"

namespace Spf {
	template<typename EngineParametersType,typename ModelType,typename RandomNumberGeneratorType>
	class AlgorithmDiag {
	public:	
		typedef typename EngineParametersType::FieldType FieldType;
		typedef std::complex<FieldType> ComplexType;
		typedef psimag::Matrix<ComplexType> MatrixType;

		AlgorithmDiag(const EngineParametersType& engineParams,ModelType& model)
			: engineParams_(engineParams),model_(model),rng_(),
					eigNew_(model.hilbertSize()),eigOld_(model.hilbertSize()),
					hilbertSize_(model_.hilbertSize()),
					matrixNew_(hilbertSize_,hilbertSize_),
					matrixOld_(hilbertSize_,hilbertSize_)
		{
		}

		void init()
		{
			model_.createHamiltonian(matrixOld_,ModelType::OLDFIELDS);
			utils::diag(matrixOld_,eigOld_,'N');
			sort(eigOld_.begin(), eigOld_.end(), std::less<FieldType>());
		}

		bool isAccepted(size_t i)
		{
			FieldType dsDirect = model_.deltaDirect(i);
				
			//FieldType oldmu=engineParams_.mu;
			
			model_.createHamiltonian(matrixNew_,ModelType::NEWFIELDS);
			diagonalize(matrixNew_,eigNew_,'N');
			//testMatrix();
			//testEigs();
			
			if (engineParams_.carriers>0) model_.adjustChemPot(eigNew_); //changes engineParams_.mu
			FieldType integrationMeasure = model_.integrationMeasure(i);
				
			return doMetropolis(dsDirect,integrationMeasure);
		}

		void accept(size_t i)
		{
			model_.accept(i);
			eigOld_ = eigNew_;
		}

		ComplexType greenFunction(size_t lambda1,size_t lambda2)
		{
			ComplexType sum = 0;
			FieldType beta = engineParams_.beta;
			FieldType mu = engineParams_.mu;
			
			for (size_t lambda=0;lambda<hilbertSize_;lambda++) 
				sum += conj(matrixNew_(lambda1,lambda)) * matrixNew_(lambda2,lambda) *utils::fermi(-beta*(eigNew_[lambda]-mu));
			//FieldType x = 0.0;
			//if (lambda1==lambda2) x = 1.0;
			return sum;
		}


		void prepare()
		{
			diagonalize(matrixNew_,eigNew_,'V',ModelType::OLDFIELDS);
		}
		
		ComplexType matrix(size_t lambda1,size_t lambda2) const
		{
			return matrixNew_(lambda1,lambda2);
		}
		
		FieldType e(size_t i) const
		{
			return eigNew_[i];
		}

		void diagonalize(MatrixType& matrix,std::vector<FieldType>& eigs,char jobz='N',size_t fields=ModelType::NEWFIELDS)
		{
			model_.createHamiltonian(matrix,fields);
			utils::diag(matrix,eigs,jobz);
			if (jobz!='V') sort(eigs.begin(), eigs.end(), std::less<FieldType>());
		}

		template<typename EngineParametersType2,typename ModelType2,typename RandomNumberGeneratorType2>
		friend std::ostream& operator<<
			(std::ostream& os,AlgorithmDiag<EngineParametersType2,ModelType2,RandomNumberGeneratorType2>& a);

	private:
		bool doMetropolis(FieldType dsDirect,FieldType integrationMeasure)
		{
			FieldType mu=engineParams_.mu;
			FieldType beta = engineParams_.beta;
			FieldType X =1.0;
			
			for (size_t i=0;i<eigNew_.size();i++) {
				FieldType temp = 0;
				if (eigNew_[i]>mu)
					temp = (1.0+exp(-beta*(eigNew_[i]-mu)))/(1.0+exp(-beta*(eigOld_[i]-mu)));
				else
				temp =(1.0+exp(beta*(eigNew_[i]-mu)))/
							(exp(beta*(eigNew_[i]-mu))+exp(-beta*(eigOld_[i]-eigNew_[i])));
			
				X *= temp;
			}
			//std::cerr<<"Xbefore="<<X<<" ";
			//if (ether.isSet("sineupdate")) X *= integrationMeasure;
			X *=  exp(-beta*dsDirect);
			X = X/(1.0+X);

			FieldType r=rng_();
			//std::cerr<<"dsDirect="<<dsDirect<<" beta="<<beta<<" X = "<<X<<" r="<<r<<"\n";
			if (X>r) return true;
			else return false;
		}
		
		void testEigs() const
		{
			FieldType eps = 1e-6;
			for (size_t i=0;i<eigOld_.size();i++) {
				if (fabs(eigOld_[i]-eigNew_[i])>eps) return;
			}
			throw std::runtime_error("Eigs are equal!!\n");
		}
		
		void testMatrix() const
		{
			FieldType eps = 1e-6;
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
		RandomNumberGeneratorType rng_;
		std::vector<FieldType> eigNew_,eigOld_;
		size_t hilbertSize_;
		MatrixType matrixNew_,matrixOld_;
		
		
	}; // AlgorithmDiag
	
	template<typename EngineParametersType,typename ModelType,typename RandomNumberGeneratorType>
	std::ostream& operator<<(std::ostream& os,AlgorithmDiag<EngineParametersType,ModelType,RandomNumberGeneratorType>& a)
	{
		
		typedef typename EngineParametersType::FieldType FieldType;
		std::vector<FieldType> eigNew(a.hilbertSize_);
		psimag::Matrix<std::complex<FieldType> > matrix(a.hilbertSize_,a.hilbertSize_);
		a.diagonalize(matrix,eigNew,'V',ModelType::OLDFIELDS);
		os<<"Eigenvalues\n";
		os<<eigNew;
		return os;
	}
} // namespace Spf

/*@}*/
#endif // ALGORITHM_DIAG_H
