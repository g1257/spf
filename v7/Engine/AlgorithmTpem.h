
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
		typedef typename PsimagLite::Vector<RealType>::Type VectorType;
		typedef MetropolisOrGlauber<RealType,RngType> MetropolisOrGlauberType;
		typedef typename EngineParametersType::IoInType IoInType;
		typedef Tpem::TpemParameters<IoInType,RealType,ModelType> TpemParametersType;
		// includes from Tpem.h
		typedef typename ModelType::MatrixType::value_type MatrixValueType;
		typedef typename ModelType::ConcurrencyType ConcurrencyType;
		typedef typename ConcurrencyType::CommType CommType;
		typedef Tpem::Tpem<TpemParametersType,MatrixValueType,ConcurrencyType> TpemType;
		typedef typename TpemType::TpemSparseType TpemSparseType;
		//
		typedef typename TpemType::ActionFunctorType ActionFunctorType;

		enum {TMPVALUES_SET,TMPVALUES_RETRIEVE};

		static const size_t computeError_ = 0;

		AlgorithmTpem(const EngineParametersType& engineParams,
		              ModelType& model,
		              IoInType& io,
					  CommType comm)
		: engineParams_(engineParams),
		  model_(model),
		  hilbertSize_(model_.hilbertSize()),
		  metropolisOrGlauber_(engineParams.detailedBalance),
		  adjustments_(engineParams),
		  adjustTpemBounds_(false),
		  tpemParameters_(io,engineParams.mu,engineParams.beta,&model),
		  tpem_(tpemParameters_,model.concurrency(),comm),
		  actionFunc_(tpemParameters_),
		  matrixOld_(model.hilbertSize(),model.hilbertSize()),
		  curMoments_(tpemParameters_.cutoff),
		  newMoments_(tpemParameters_.cutoff),
		  error_("DISABLED")
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
			newMoments_ = curMoments_ + moments;
			
			adjustChemPot(newMoments_,i,engineParams_.adjustEach);
			//RealType integrationMeasure = model_.integrationMeasure(i);
			RealType X = exp(dS);
			return metropolisOrGlauber_(X,rng);
		}

		void accept(size_t i)
		{
			model_.accept(i);
			// update current moments
			curMoments_ = newMoments_;
			matrixOld_ = matrixNew_;
		}

		void prepare()
		{
			TpemSparseType matrix(hilbertSize_,hilbertSize_);
			setMatrix(matrix,ModelType::OLDFIELDS);
			assert(isHermitian(matrix));
			tpem_.calcMoments(matrix,curMoments_);
			adjustChemPot(curMoments_,0,0);
// 			std::cerr<<"IsHerm = "<<b<<"\n";
		}
		
		std::string error() const
		{
			return error_;
		}

		const VectorType& moment() const { return curMoments_; }

		ModelType& model() { return model_; }
		
		TpemType& tpem() { return tpem_; }

		const TpemParametersType& tpemParameters() const { return tpemParameters_; }

		template<typename EngineParametersType2,typename ModelType2,
		         typename RandomNumberGeneratorType2>
		friend std::ostream& operator<<(std::ostream& os,
                                        const AlgorithmTpem<EngineParametersType2,
		                                ModelType2,RandomNumberGeneratorType2>& a);

	private:

		void adjustChemPot(const VectorType& moments,size_t i,size_t adjustEach)
		{
			if (engineParams_.carriers<=0) return;
			if (!needsAdjustment(i,adjustEach)) return;

			try {
				engineParams_.mu = adjustments_.adjChemPot(moments,tpem_);
			} catch (std::exception& e) {
				std::cerr<<e.what()<<"\n";
			}
		}

		bool needsAdjustment(size_t i,size_t adjustEach) const
		{
			if (adjustEach==0) return true;
			size_t x = (i % adjustEach);
			return (x==0);
		}

		RealType calcDeltaAction(VectorType& moments,
		                         const TpemSparseType& matrix0,
		                         const TpemSparseType& matrix1)
		{
			VectorType actionCoeffs(tpemParameters_.cutoff);
			tpem_.calcCoeffs(actionCoeffs,actionFunc_); 

			tpem_.calcMomentsDiff(moments,matrix0, matrix1);


			size_t total = moments.size();
			size_t start = size_t(total*0.8);
			
			if (!computeError_ || start<2)  {
				return -tpem_.expand(actionCoeffs, moments);
			}
		
			RealType error1 = 0;
			RealType x = -tpem_.expand(actionCoeffs, moments,0,start-1);
			for (size_t i=start;i<total;i++) {
				RealType deltaX = -tpem_.expand(actionCoeffs, moments,i-1,i);
				error1 += deltaX * deltaX;
				x += deltaX;
			}
			error_ = ttos(error1);
			//assert(x==-tpem_.expand(actionCoeffs, moments));
			
			return x;
		}

		void setMatrix(TpemSparseType& matrix,size_t oldOrNewFields) const
		{
			model_.createHsparse(matrix,oldOrNewFields);
			TpemSparseType diagB(matrix.row(),matrix.col());
			diagB.makeDiagonal(matrix.row(),-tpemParameters_.b);
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
		std::string error_;
	}; // AlgorithmTpem
	
	template<typename EngineParametersType,typename ModelType,
	         typename RandomNumberGeneratorType>
	std::ostream& operator<<(std::ostream& os,const AlgorithmTpem<
	    EngineParametersType,ModelType,RandomNumberGeneratorType>& a)
	{
		os<<a.tpemParameters_;
		return os;
	}
} // namespace Spf

/*@}*/
#endif // ALGORITHM_TPEM_H
