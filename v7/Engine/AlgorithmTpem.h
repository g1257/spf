
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

namespace Spf {
	template<typename EngineParametersType,typename ModelType,typename RngType>
	class AlgorithmTpem {
	public:	
		typedef typename EngineParametersType::FieldType RealType;
		typedef std::complex<RealType> ComplexType;
		typedef PsimagLite::Matrix<ComplexType> MatrixType;
		typedef std::vector<RealType> VectorType;
		typedef tpem_sparse TpemSparseType;
		typedef TpemOptions TpemOptionsType;

		enum {TMPVALUES_SET,TMPVALUE_RETRIEVE};

		AlgorithmTpem(const EngineParametersType& engineParams,ModelType& model)
		: engineParams_(engineParams),model_(model),
		  hilbertSize_(model_.hilbertSize()),
		  actionCoeffs_(cutoff_),
		  tpemOptions_(),
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
			if (engineParams_.isSet("adjusttpembounds")) {
				if (tpemAdjustBounds(aux.sparseTmp[0],et)!=0) {
					cerr<<"Cannot adjust bounds for tpem spectrum\n";
					exit(1);
				}
			}
			RealType dS = calcDeltaAction(moment_,aux.sparseMatrix[0], aux.sparseTmp[0],
			                 support);
			
			dS -= engineParams_.beta*dsDirect;

			
			//if (engineParams_.carriers>0) model_.adjustChemPot(eigNew_); //changes engineParams_.mu
			//RealType integrationMeasure = model_.integrationMeasure(i);
				
			return metropolisOrGlauber(dsDirect,rng);
		}

		void accept(size_t i)
		{
			model_.accept(i);
			// update current moments
			for (size_t j=0;j<cutoff_;j++) curMoments_[j] -= moment[j];
		}

		void prepare()
		{
// 			diagonalize(matrixNew_,eigNew_,'V',ModelType::OLDFIELDS);
		}
		
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
		
		double calcDeltaAction(const VectorType &moment,
		                       const TpemSparseType* matrix1,
		                       const TpemSparseType* matrix, 
		                       const std::vector<size_t>& support)
		{
			
			RealType e1=aux.varTpem_b-aux.varTpem_a;
			RealType e2=aux.varTpem_a+aux.varTpem_b;
			RealType e11=aux.btmp-aux.atmp;
			RealType e21=aux.atmp+aux.btmp;

			if (e1 > e11) e1 = e11;
			if (e2 < e21) e2 = e21;

			if (engineParams_.carriers>0) {
				tmpValues(a,b,mu,beta,0);
				tpem_calculate_coeffs (actionCoeffs_,actionFunc_,tpemOptions_);
			} 
			if (engineParams_.isSet("adjusttpembounds")) {
				a=0.5*(e2-e1);
				b=0.5*(e2+e1);
				tpem_calculate_coeffs (actionCoeffs_,actionFunc_,tpemOptions_);
			}
		
			tmpValues(a,b,mu,beta,0);
			
			TpemSparseType* mod_matrix1 = matrix1;
			TpemSparseType* mod_matrix = matrix;

			if (engineParams_.isSet("adjusttpembounds")) {
					// full copy/allocation
					mod_matrix1 = new_tpem_sparse(engineParams_.hilbertSize,matrix1->rowptr[engineParams_.hilbertSize]);
					tpem_sparse_copy (matrix1, mod_matrix1);
					mod_matrix = new_tpem_sparse(engineParams_.hilbertSize,matrix->rowptr[engineParams_.hilbertSize]);
					tpem_sparse_copy (matrix, mod_matrix);
					
					tpem_sparse_scale(mod_matrix1,a/aux.atmp,(b-aux.btmp)*aux.atmp/(a*a));
					tpem_sparse_scale(mod_matrix,a/aux.varTpem_a,(b-aux.varTpem_b)*aux.varTpem_a/(a*a));
					
			}

			tpem_calculate_moment_diff (mod_matrix1, mod_matrix,moment_,support, tpemOptions_);

			dS = -tpem_expansion (coeffs, moment_);

			if (engineParams_.isSet("adjusttpembounds")) {
				tpem_sparse_free(mod_matrix);
				tpem_sparse_free(mod_matrix1);
			}
			return dS;
		}

		void metropolisOrGlauber(const RealType& dS) const
		{
			if (DO_GLAUBER) {
				if (dS<0) {
					dS=exp(dS)/(1.0+exp(dS));
				} else {
					dS=1.0/(1.0+exp(-dS));
				}
				return (dS>myRandom());
			}
			// METROPOLIS PROPER 
			return (dS > 0.0 || myRandom () < exp (dS));
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
		size_t cutoff_;
		RealType a_,b_,mu_,beta_;
		TpemOptionsType tpemOptions_;
		VectorType actionCoeffs_;
		VectorType moment0_;
		VectorType moment1_;
	}; // AlgorithmTpem
	
	template<typename EngineParametersType,typename ModelType,
		typename RandomNumberGeneratorType>
	std::ostream& operator<<(std::ostream& os,const AlgorithmTpem<
			EngineParametersType,ModelType,RandomNumberGeneratorType>& a)
	{
		
		typedef typename EngineParametersType::RealType RealType;
		std::vector<RealType> eigNew(a.hilbertSize_);
		PsimagLite::Matrix<std::complex<RealType> > matrix(a.hilbertSize_,a.hilbertSize_);
		a.diagonalize(matrix,eigNew,'V',ModelType::OLDFIELDS);
		os<<"Eigenvalues\n";
		os<<eigNew;
		return os;
	}
} // namespace Spf

/*@}*/
#endif // ALGORITHM_TPEM_H
