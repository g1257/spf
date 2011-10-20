
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
#include "Adjustments.h"

namespace Spf {
	template<typename EngineParametersType,typename ModelType,typename RngType>
	class AlgorithmTpem {
		static const bool DO_GLAUBER = true;
		typedef Adjustments<EngineParametersType> AdjustmentsType;
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
		  adjustments_(engineParams),
		  adjustTpemBounds_(false),
		  tpemParameters_(io,engineParams.mu,engineParams.beta,&model),
		  tpem_(tpemParameters_),
		  actionFunc_(tpemParameters_),
		  matrixOld_(model.hilbertSize(),model.hilbertSize()),
		  curMoments_(tpemParameters_.cutoff),
		  newMoments_(tpemParameters_.cutoff)
		{}

		void init()
		{
 			setMatrix(matrixOld_,ModelType::OLDFIELDS);
			tpem_.calcMoments(matrixOld_,curMoments_);
		}

		size_t hilbertSize() const { return hilbertSize_; }

		bool isAccepted(size_t i,RngType& rng)
		{
			RealType dsDirect = model_.deltaDirect(i);
							
			setMatrix(matrixNew_,ModelType::NEWFIELDS);

			VectorType moments(curMoments_.size());
			RealType dS = calcDeltaAction(moments,matrixOld_, matrixNew_);

			dS -= engineParams_.beta*dsDirect;
			newMoments_ = curMoments_ - moments;
			
			adjustChemPot(newMoments_);
			//RealType integrationMeasure = model_.integrationMeasure(i);
			RealType X = exp(dS);
			return metropolisOrGlauber_(X,rng);
		}

		void accept(size_t i)
		{
			model_.accept(i);
			// update current moments
			curMoments_ = newMoments_;
		}

		void prepare()
		{
			TpemSparseType matrix(hilbertSize_,hilbertSize_);
			setMatrix(matrix,ModelType::OLDFIELDS);
			assert(isHermitian(matrix));
			tpem_.calcMoments(matrix,curMoments_);
			adjustChemPot(curMoments_);
// 			std::cerr<<"IsHerm = "<<b<<"\n";
		}
		
		void adjustChemPot(const VectorType& moments)
		{
			if (engineParams_.carriers<=0) return;
			try {
				engineParams_.mu = adjustments_.adjChemPot(moments,tpem_); //changes engineParams_.mu
			} catch (std::exception& e) {
				std::cerr<<e.what()<<"\n";
			}
		}
		
		const VectorType& moment() const { return curMoments_; }

		ModelType& model() { return model_; }
		
		TpemType& tpem() { return tpem_; }

		const TpemParametersType& tpemParameters() const { return tpemParameters_; }

		template<typename EngineParametersType2,typename ModelType2,
			typename RandomNumberGeneratorType2>
		friend std::ostream& operator<<
			(std::ostream& os,const AlgorithmTpem<EngineParametersType2,
					ModelType2,RandomNumberGeneratorType2>& a);

	private:

		double calcDeltaAction(VectorType& moments,
		                       const TpemSparseType& matrix0,
		                       const TpemSparseType& matrix1)
		{
			VectorType actionCoeffs(tpemParameters_.cutoff);
			tpem_.calcCoeffs(actionCoeffs,actionFunc_); 
			
			tpem_.calcMomentsDiff(moments,matrix0, matrix1);

			return -tpem_.expand(actionCoeffs, moments);
		}

		void setMatrix(TpemSparseType& matrix,size_t oldOrNewFields) const
		{
			model_.createHsparse(matrix,oldOrNewFields);
			TpemSparseType diagB(matrix.rank(),matrix.rank());
			diagB.makeDiagonal(matrix.rank(),-tpemParameters_.b);
			matrix += diagB;
			matrix *= (1.0/tpemParameters_.a);
		}

		const EngineParametersType& engineParams_;
		ModelType& model_;
		size_t hilbertSize_;
		MetropolisOrGlauberType metropolisOrGlauber_;
		AdjustmentsType adjustments_;
		bool adjustTpemBounds_;
		TpemParametersType tpemParameters_;
		TpemType tpem_;
		ActionFunctorType actionFunc_;
		TpemSparseType matrixOld_;
		TpemSparseType matrixNew_;
		VectorType curMoments_;
		VectorType newMoments_;
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
