
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
#include "Tpem.h"
#include "MetropolisOrGlauber.h"
#include "TpemParameters.h"

namespace Spf {
	template<typename EngineParametersType,typename ModelType,typename RngType>
	class AlgorithmTpem {
		static const bool DO_GLAUBER = true;

	public:	
		typedef typename EngineParametersType::RealType RealType;
		typedef std::complex<RealType> ComplexType;
		typedef PsimagLite::Matrix<ComplexType> MatrixType;
		typedef std::vector<RealType> VectorType;
		typedef MetropolisOrGlauber<RealType,RngType> MetropolisOrGlauberType;
		typedef typename EngineParametersType::IoInType IoInType;
		typedef Tpem::TpemParameters<IoInType,RealType,ModelType> TpemParametersType;
		// includes from Tpem.h
		typedef Tpem::Tpem<TpemParametersType,typename ModelType::MatrixType::value_type> TpemType;
		typedef typename TpemType::TpemSparseType TpemSparseType;
		//
		typedef typename TpemType::ActionFunctorType ActionFunctorType;
		
		
		enum {TMPVALUES_SET,TMPVALUES_RETRIEVE};

		AlgorithmTpem(const EngineParametersType& engineParams,
		              ModelType& model,
		              IoInType& io)
		: engineParams_(engineParams),
		  model_(model),
		  hilbertSize_(model_.hilbertSize()),
		  metropolisOrGlauber_(),
		  adjustTpemBounds_(false),
		  tpemParameters_(io,engineParams.mu,engineParams.beta,&model),
		  tpem_(tpemParameters_),
		  actionFunc_(tpemParameters_),
		  matrixOld_(model.hilbertSize(),model.hilbertSize()),
		  actionCoeffs_(tpemParameters_.cutoff),
		  moment_(tpemParameters_.cutoff),
		  curMoments_(tpemParameters_.cutoff)
		{
			tpem_.calcCoeffs(actionCoeffs_,actionFunc_); 
		}

		void init()
		{
 			setMatrix(matrixOld_,ModelType::OLDFIELDS);
// 			diag(matrixOld_,eigOld_,'N');
// 			sort(eigOld_.begin(), eigOld_.end(), std::less<RealType>());
		}

		size_t hilbertSize() const { return hilbertSize_; }

		bool isAccepted(size_t i,RngType& rng)
		{
			RealType dsDirect = model_.deltaDirect(i);
							
			setMatrix(matrixNew_,ModelType::NEWFIELDS);

			RealType dS = calcDeltaAction(matrixOld_, matrixNew_);
			
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
			for (size_t j=0;j<tpemParameters_.cutoff;j++)
				curMoments_[j] -= moment_[j];
		}

		void prepare()
		{
			TpemSparseType matrix(hilbertSize_,hilbertSize_);
			setMatrix(matrix,ModelType::OLDFIELDS);
			bool b = isHermitian(matrix,true);
			tpem_.calcMoments(matrix,moment_);
			std::cerr<<"IsHerm = "<<b<<"\n";
		}
		
		const VectorType& moment() const { return moment_; }

		ModelType& model() { return model_; }
		
		TpemType& tpem() { return tpem_; }

		const TpemParametersType& tpemParameters() const { return tpemParameters_; }

		template<typename EngineParametersType2,typename ModelType2,
			typename RandomNumberGeneratorType2>
		friend std::ostream& operator<<
			(std::ostream& os,const AlgorithmTpem<EngineParametersType2,
					ModelType2,RandomNumberGeneratorType2>& a);

	private:
		
		double calcDeltaAction(TpemSparseType& matrix1,
		                       TpemSparseType& matrix)
		{
			tpem_.calcMomentsDiff(moment_,matrix1, matrix);

			return -tpem_.expand(actionCoeffs_, moment_);
		}

		void setMatrix(TpemSparseType& matrix,size_t oldOrNewFields) const
		{
			model_.createHsparse(matrix,oldOrNewFields);
			
			typedef typename TpemSparseType::value_type ElementType;
			size_t counter = 0;
			for (size_t i=0;i<matrix.rank();i++) {
				for (int k=matrix.getRowPtr(i);k<matrix.getRowPtr(i+1);k++) {
					counter++;
					if (size_t(matrix.getCol(k))!=i) continue;
					ElementType tmp = matrix.getValue(k)-tpemParameters_.b;
					matrix.setValues(counter-1,tmp);
				}
			}
			matrix *= (1.0/tpemParameters_.a);
		}

		const EngineParametersType& engineParams_;
		ModelType& model_;
		size_t hilbertSize_;
		MetropolisOrGlauberType metropolisOrGlauber_;
		bool adjustTpemBounds_;
		TpemParametersType tpemParameters_;
		TpemType tpem_;
		ActionFunctorType actionFunc_;
		TpemSparseType matrixOld_;
		TpemSparseType matrixNew_;
		VectorType actionCoeffs_;
		VectorType moment_;
		VectorType curMoments_;
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
		os<<a.tpemParameters_;
		os<<"operator<< (os,AlgorithmTpem) unimplemented yet(sorry)\n";
		return os;
	}
} // namespace Spf

/*@}*/
#endif // ALGORITHM_TPEM_H
