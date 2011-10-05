
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

namespace Spf {
	template<typename EngineParametersType,typename ModelType,typename RngType>
	class AlgorithmTpem {
	public:	
		typedef typename EngineParametersType::FieldType FieldType;
		typedef std::complex<FieldType> ComplexType;
		typedef PsimagLite::Matrix<ComplexType> MatrixType;

		AlgorithmTpem(const EngineParametersType& engineParams,ModelType& model)
		: engineParams_(engineParams),model_(model),
		  hilbertSize_(model_.hilbertSize()),
		  actionCoeffs_(cutoff_),
		  moment0_(cutoff_),
		  moment1_(cutoff_)
		{
			tmpValues(a,b,mu,beta,0);
			tpem_calculate_coeffs (actionCoeffs_,actionFunc_,tpemOptions_); 
		}

		void init()
		{
// 			model_.createHamiltonian(matrixOld_,ModelType::OLDFIELDS);
// 			diag(matrixOld_,eigOld_,'N');
// 			sort(eigOld_.begin(), eigOld_.end(), std::less<FieldType>());
		}

		size_t hilbertSize() const { return hilbertSize_; }

		bool isAccepted(size_t i,RngType& rng)
		{
			FieldType dsDirect = model_.deltaDirect(i);
							
			model_.createHsparse(matrixNew_,ModelType::NEWFIELDS);
			if (ether.isSet("adjusttpembounds")) {
				if (tpemAdjustBounds(aux.sparseTmp[0],ether,aux)!=0) {
					cerr<<"Cannot adjust bounds for tpem spectrum\n";
					exit(1);
				}
			}
			dS += calcDeltaAction(moment,aux.sparseMatrix[0], aux.sparseTmp[0],
								support,ether,aux,tpemOptions);
			
			dS -= ether.beta*dsDirect;

			
			//if (engineParams_.carriers>0) model_.adjustChemPot(eigNew_); //changes engineParams_.mu
			//FieldType integrationMeasure = model_.integrationMeasure(i);
				
			return metropolisOrGlauber(dsDirect,rng);
		}

		void accept(size_t i)
		{
			model_.accept(i);
			// update current moments
			for (j=0;j<ether.tpem_cutoff;j++) aux.curMoments[j] -= moment[j];
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
// 		const FieldType& e(size_t i) const
// 		{
// 			return eigNew_[i];
// 		}
// 
// 		void diagonalize(
// 				MatrixType& matrix,
// 				std::vector<FieldType>& eigs,
// 				char jobz='N',
// 				size_t fields=ModelType::NEWFIELDS) const
// 		{
// 			model_.createHamiltonian(matrix,fields);
// 			diag(matrix,eigs,jobz);
// 			if (jobz!='V')
// 				std::sort(eigs.begin(), eigs.end(), std::less<FieldType>());
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
		
		double calcDeltaAction(vector<double> &moment,tpem_sparse *matrix1, tpem_sparse *matrix, 
							 vector<unsigned int> const &support,Parameters const &ether,Aux &aux,TpemOptions const &tpemOptions)
		{
			
			FieldType e1=aux.varTpem_b-aux.varTpem_a;
			FieldType e2=aux.varTpem_a+aux.varTpem_b;
			FieldType e11=aux.btmp-aux.atmp;
			FieldType e21=aux.atmp+aux.btmp;
			
			if (e1 > e11) e1 = e11;
			if (e2 < e21) e2 = e21;
			
			
			if (ether.carriers>0) {
				tmpValues(a,b,mu,beta,0);
				tpem_calculate_coeffs (actionCoeffs_,actionFunc_,tpemOptions_);
			} 
			if (ether.isSet("adjusttpembounds")) {
				a=0.5*(e2-e1);
				b=0.5*(e2+e1);
				tpem_calculate_coeffs (actionCoeffs_,actionFunc_,tpemOptions_);
			}
		
			tmpValues(a,b,mu,beta,0);
			
			tpem_sparse *mod_matrix1, *mod_matrix;
			
			if (ether.isSet("adjusttpembounds")) {
					// full copy/allocation
					mod_matrix1 = new_tpem_sparse(ether.hilbertSize,matrix1->rowptr[ether.hilbertSize]);
					tpem_sparse_copy (matrix1, mod_matrix1);
					mod_matrix = new_tpem_sparse(ether.hilbertSize,matrix->rowptr[ether.hilbertSize]);
					tpem_sparse_copy (matrix, mod_matrix);
					
					tpem_sparse_scale(mod_matrix1,a/aux.atmp,(b-aux.btmp)*aux.atmp/(a*a));
					tpem_sparse_scale(mod_matrix,a/aux.varTpem_a,(b-aux.varTpem_b)*aux.varTpem_a/(a*a));
					
			} else {		
				// just pointers
				mod_matrix1 = matrix1;
				mod_matrix = matrix;
			}
			
			tpem_calculate_moment_diff (mod_matrix1, mod_matrix,moment,support, tpemOptions);
			
			dS = -tpem_expansion (coeffs, moment);
			
			if (ether.isSet("adjusttpembounds")) {
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

		const EngineParametersType& engineParams_;
		ModelType& model_;
		size_t hilbertSize_;
	}; // AlgorithmTpem
	
	template<typename EngineParametersType,typename ModelType,
		typename RandomNumberGeneratorType>
	std::ostream& operator<<(std::ostream& os,const AlgorithmTpem<
			EngineParametersType,ModelType,RandomNumberGeneratorType>& a)
	{
		
		typedef typename EngineParametersType::FieldType FieldType;
		std::vector<FieldType> eigNew(a.hilbertSize_);
		PsimagLite::Matrix<std::complex<FieldType> > matrix(a.hilbertSize_,a.hilbertSize_);
		a.diagonalize(matrix,eigNew,'V',ModelType::OLDFIELDS);
		os<<"Eigenvalues\n";
		os<<eigNew;
		return os;
	}
} // namespace Spf

/*@}*/
#endif // ALGORITHM_TPEM_H
