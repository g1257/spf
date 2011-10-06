
/** \ingroup SPF */
/*@{*/

/*! \file AlgorithmTpem.h
 *
 *  Diagonalization method for SPF
 *
 */
#ifndef ALGORITHM_TPEM_H
#define ALGORITHM_TPEM_H
#include <algorithm>
#include "ProgressIndicator.h" // in PsimagLite
#include "Matrix.h" // in PsimagLite
#include "Fermi.h" // in PsimagLite
#include "Complex.h" // in PsimagLite
#include "tpemplus.h"
#include "TpemFunctors.h"
#include "MetropolisOrGlauber.h"

namespace Spf {
	template<typename EngineParametersType,typename ModelType,typename RngType>
	class AlgorithmTpem {
		static const bool DO_GLAUBER = true;

	public:	
		typedef typename EngineParametersType::FieldType RealType;
		typedef std::complex<RealType> ComplexType;
		typedef PsimagLite::Matrix<ComplexType> MatrixType;
		typedef std::vector<RealType> VectorType;
		typedef MetropolisOrGlauber<RealType,RngType> MetropolisOrGlauberType;
		// includes from tpemplus.h:
		typedef tpem_sparse TpemSparseType;
		typedef TpemOptions TpemOptionsType;
		// 
		typedef ObservableFunctor ObservableFunctorType;
		
		enum {TMPVALUES_SET,TMPVALUE_RETRIEVE};

		AlgorithmTpem(const EngineParametersType& engineParams,ModelType& model)
		: engineParams_(engineParams),model_(model),
		  hilbertSize_(model_.hilbertSize()),
		  metropolisOrGlauber_(),
		  tpemOptions_(),
		  adjustTpemBounds_(false),
		  actionCoeffs_(cutoff_),
		  moment0_(cutoff_),
		  moment1_(cutoff_)
		{
			throw std::runtime_error("Need to set cutoff_, a_, b_, mu_, beta_\n");
			//tmpValues(a,b,mu,beta,TMPVALUES_SET);
			tpem_calculate_coeffs (actionCoeffs_,actionFunc_,tpemOptions_); 
		}

		void init()
		{
// 			model_.createHamiltonian(matrixOld_,ModelType::OLDFIELDS);
// 			diag(matrixOld_,eigOld_,'N');
// 			sort(eigOld_.begin(), eigOld_.end(), std::less<RealType>());
		}

		size_t hilbertSize() const { return hilbertSize_; }

		bool isAccepted(size_t i,RngType& rng)
		{
			RealType dsDirect = model_.deltaDirect(i);
							
			model_.createHsparse(matrixNew_,ModelType::NEWFIELDS);
			if (adjustTpemBounds_) {
				if (tpemAdjustBounds(matrixNew_)!=0) {
					throw std::runtime_error(
						"Cannot adjust bounds for tpem spectrum\n");
				}
			}
			RealType dS = calcDeltaAction(moment_,matrixOld_, matrixNew_);
			
			dS -= engineParams_.beta*dsDirect;

			
			//if (engineParams_.carriers>0) model_.adjustChemPot(eigNew_); //changes engineParams_.mu
			//RealType integrationMeasure = model_.integrationMeasure(i);
			RealType X = exp(dS);
			return metropolisOrGlauber_(X,rng);
		}

		void accept(size_t i)
		{
			model_.accept(i);
			// update current moments
			for (size_t j=0;j<cutoff_;j++) curMoments_[j] -= moment_[j];
		}

		void prepare()
		{
// 			diagonalize(matrixNew_,eigNew_,'V',ModelType::OLDFIELDS);
		}
		
		const VectorType& moment() const { return moment_; }
		
		const TpemOptions& tpemOptions() const { return tpemOptions_; }
		
		const size_t cutoff() const { return cutoff_; }

		ModelType& model() { return model_; }

// 		const ComplexType& matrix(size_t lambda1,size_t lambda2) const
// 		{
// 			return matrixNew_(lambda1,lambda2);
// 		}
// 		
// 		const RealType& e(size_t i) const
// 		{
// 			return eigNew_[i];
// 		}
// 
// 		void diagonalize(
// 				MatrixType& matrix,
// 				std::vector<RealType>& eigs,
// 				char jobz='N',
// 				size_t fields=ModelType::NEWFIELDS) const
// 		{
// 			model_.createHamiltonian(matrix,fields);
// 			diag(matrix,eigs,jobz);
// 			if (jobz!='V')
// 				std::sort(eigs.begin(), eigs.end(), std::less<RealType>());
// 		}
// 		
// 		void printMatrix(size_t mode) const
// 		{
// 			if (mode==ModelType::NEWFIELDS) {
// 				std::cerr<<matrixNew_;
// 				return;
// 			}
// 			MatrixType m(matrixNew_.n_row(),matrixNew_.n_col());
// 			if (!isHermitian(m)) throw std::runtime_error("Problem\n");
// 			model_.createHamiltonian(m,ModelType::OLDFIELDS);
// 			std::cerr<<m;
// 		}

		template<typename EngineParametersType2,typename ModelType2,
			typename RandomNumberGeneratorType2>
		friend std::ostream& operator<<
			(std::ostream& os,const AlgorithmTpem<EngineParametersType2,
					ModelType2,RandomNumberGeneratorType2>& a);

	private:
		
		int tpemAdjustBounds(const TpemSparseType* matrix) const
		{
			// unimplemented for now
			return 0;
		}
		
		double calcDeltaAction(const VectorType &moment,
		                       const TpemSparseType* matrix1,
		                       const TpemSparseType* matrix)
		{
			RealType a = a_;
			RealType b = b_;
			RealType mu = mu_;
			RealType beta = beta_;
			RealType e1=b_-a_;
			RealType e2=b_+a_;
			RealType e11=engineParams_.b - engineParams_.b;
			RealType e21=engineParams_.b + engineParams_.b;

			if (e1 > e11) e1 = e11;
			if (e2 < e21) e2 = e21;

			if (engineParams_.carriers>0) {
				tmpValues(a,b,mu,beta,TMPVALUES_SET);
				tpem_calculate_coeffs (actionCoeffs_,actionFunc_,tpemOptions_);
			} 
			if (adjustTpemBounds_) {
				a=0.5*(e2-e1);
				b=0.5*(e2+e1);
				tpem_calculate_coeffs (actionCoeffs_,actionFunc_,tpemOptions_);
			}
		
			tmpValues(a,b,mu,beta,TMPVALUES_SET);
			
			TpemSparseType* mod_matrix1 = matrix1;
			TpemSparseType* mod_matrix = matrix;

			if (adjustTpemBounds_) {
				std::string s(__FILE__);
				s += " " + __FUNCTION__ + " adjusttpembounds, not sure how to proceed\n";
				throw std::runtime_error(s.c_str());
					// full copy/allocation
// 					mod_matrix1 = new_tpem_sparse(engineParams_.hilbertSize,matrix1->rowptr[engineParams_.hilbertSize]);
// 					tpem_sparse_copy (matrix1, mod_matrix1);
// 					mod_matrix = new_tpem_sparse(engineParams_.hilbertSize,matrix->rowptr[engineParams_.hilbertSize]);
// 					tpem_sparse_copy (matrix, mod_matrix);
// 					
// 					tpem_sparse_scale(mod_matrix1,a/aux.atmp,(b-aux.btmp)*aux.atmp/(a*a));
// 					tpem_sparse_scale(mod_matrix,a/aux.varTpem_a,(b-aux.varTpem_b)*aux.varTpem_a/(a*a));
					
			}

			tpem_calculate_moment_diff (mod_matrix1, mod_matrix,moment_,support_, tpemOptions_);

			RealType dS = -tpem_expansion (actionCoeffs_, moment_);

			if (adjustTpemBounds_) {
				tpem_sparse_free(mod_matrix);
				tpem_sparse_free(mod_matrix1);
			}
			return dS;
		}
		
		void tmpValues(RealType &a,RealType &b,RealType &mu,RealType &beta,size_t option)
		{
			if (option==TMPVALUES_SET) { // set static members
				a_ = a;
				b_ = b;
				mu_ = mu;
				beta_ = beta;
				return;
			}
			//else retrieve static members
			a = a_;
			b = b_;
			mu = mu_;
			beta = beta_;
		}

		const EngineParametersType& engineParams_;
		ModelType& model_;
		size_t hilbertSize_;
		MetropolisOrGlauberType metropolisOrGlauber_;
		size_t cutoff_;
		RealType a_,b_,mu_,beta_;
		TpemOptionsType tpemOptions_;
		std::vector<size_t> support_;
		bool adjustTpemBounds_;
		ObservableFunctor actionFunc_;
		TpemSparseType* matrixOld_;
		TpemSparseType* matrixNew_;
		VectorType actionCoeffs_;
		VectorType curMoments_;
		VectorType moment_;
		VectorType moment0_;
		VectorType moment1_;
	}; // AlgorithmTpem
	
	template<typename EngineParametersType,typename ModelType,
		typename RandomNumberGeneratorType>
	std::ostream& operator<<(std::ostream& os,const AlgorithmTpem<
			EngineParametersType,ModelType,RandomNumberGeneratorType>& a)
	{
		
// 		typedef typename EngineParametersType::FieldType RealType;
// 		std::vector<RealType> eigNew(a.hilbertSize_);
// 		PsimagLite::Matrix<std::complex<RealType> > matrix(a.hilbertSize_,a.hilbertSize_);
// 		a.diagonalize(matrix,eigNew,'V',ModelType::OLDFIELDS);
// 		os<<"Eigenvalues\n";
// 		os<<eigNew;
		os<<"operator<< (os,AlgorithmTpem) unimplemented yet(sorry)\n";
		return os;
	}
} // namespace Spf

/*@}*/
#endif // ALGORITHM_TPEM_H
