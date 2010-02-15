
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
					matrix_(hilbertSize_,hilbertSize_),needsDiagonalization_(true)
		{
		}

		void init()
		{
			model_.createHamiltonian(matrix_,ModelType::OLDFIELDS);
			utils::diag(matrix_,eigOld_,'N');
			sort(eigOld_.begin(), eigOld_.end(), std::less<FieldType>());
		}

		bool isAccepted(size_t i)
		{
			FieldType dsDirect = model_.deltaDirect(i);
				
			//FieldType oldmu=engineParams_.mu;
			
			model_.createHamiltonian(matrix_,ModelType::NEWFIELDS);
			diagonalize(matrix_,eigNew_,'N');
			
			model_.adjustChemPot(eigNew_); //changes engineParams_.mu
			FieldType integrationMeasure = model_.integrationMeasure(i);
				
			return doMetropolis(dsDirect,integrationMeasure);
		}

		void accept(size_t i)
		{
			model_.accept(i);
			eigOld_ = eigNew_;
			needsDiagonalization_ = true;
		}

		ComplexType greenFunction(size_t lambda1,size_t lambda2)
		{
			ComplexType sum = 0;
			FieldType beta = engineParams_.beta;
			FieldType mu = engineParams_.mu;
			if (needsDiagonalization_) diagonalize(matrix_,eigNew_,'V');
			needsDiagonalization_ = false;
			return matrix_(lambda1,lambda2);
			
			for (size_t lambda=0;lambda<hilbertSize_;lambda++) 
				sum += conj(matrix_(lambda1,lambda)) * matrix_(lambda2,lambda) *utils::fermi(beta*(eigNew_[lambda]-mu));
			//FieldType x = 0.0;
			//if (lambda1==lambda2) x = 1.0;
			return sum;
		}
		
		ComplexType matrix(size_t lambda1,size_t lambda2)
		{
			ComplexType sum = 0;
			//FieldType beta = engineParams_.beta;
			//FieldType mu = engineParams_.mu;
			if (needsDiagonalization_) diagonalize(matrix_,eigNew_,'V');
			needsDiagonalization_ = false;
			return matrix_(lambda1,lambda2);
			
			/*for (size_t lambda=0;lambda<hilbertSize_;lambda++) 
				sum += conj(matrix_(lambda1,lambda)) * matrix_(lambda2,lambda) *utils::fermi(beta*(eigNew_[lambda]-mu));
			//FieldType x = 0.0;
			//if (lambda1==lambda2) x = 1.0;
			return sum;*/
		}
		
		FieldType e(size_t i)
		{
			if (needsDiagonalization_) diagonalize(matrix_,eigNew_,'V');
			needsDiagonalization_ = false;
			return eigNew_[i];
		}

		void diagonalize(MatrixType& matrix,std::vector<FieldType>& eigs,char jobz='N')
		{
			model_.createHamiltonian(matrix_,ModelType::NEWFIELDS);
			utils::diag(matrix_,eigNew_,jobz);
			if (jobz!='V') sort(eigNew_.begin(), eigNew_.end(), std::less<FieldType>());
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
					temp = (double)(1.0+exp(-beta*(eigNew_[i]-mu)))/(1.0+exp(-beta*(eigOld_[i]-mu)));
				else
				temp =(double)(1.0+exp(beta*(eigNew_[i]-mu)))/
							(exp(beta*(eigNew_[i]-mu))+exp(-beta*(eigOld_[i]-eigNew_[i])));
			
				X *= temp;
			}

			//if (ether.isSet("sineupdate")) X *= integrationMeasure;
			X *=  exp(-beta*dsDirect);
			X = X/(1.0+X);

			FieldType r=rng_();

			if (X>r) return true;
			else return false;
		}

		const EngineParametersType& engineParams_;
		ModelType& model_;
		RandomNumberGeneratorType rng_;
		std::vector<FieldType> eigNew_,eigOld_;
		size_t hilbertSize_;
		MatrixType matrix_;
		bool needsDiagonalization_;
		
		
	}; // AlgorithmDiag
	
	template<typename EngineParametersType,typename ModelType,typename RandomNumberGeneratorType>
	std::ostream& operator<<(std::ostream& os,AlgorithmDiag<EngineParametersType,ModelType,RandomNumberGeneratorType>& a)
	{
		
		typedef typename EngineParametersType::FieldType FieldType;
		std::vector<FieldType> eigNew(a.hilbertSize_);
		psimag::Matrix<std::complex<FieldType> > matrix(a.hilbertSize_,a.hilbertSize_);
		a.diagonalize(matrix,eigNew,'V');
		os<<"Eigenvalues\n";
		os<<eigNew;
		return os;
	}
} // namespace Spf

/*@}*/
#endif // ALGORITHM_DIAG_H
