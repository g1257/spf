
#include "TpemParameters.h"
#include "Tpem.h"
#include "RandomForTests.h"
#include "IoSimple.h"
#include "Concurrency.h"
#include "Vector.h"

template<typename RealType>
typename PsimagLite::EnableIf<!PsimagLite::IsComplexNumber<RealType>::True,RealType>::Type
realOrComplexHopping()
{
	return -1;
}

template<typename ComplexType>
typename PsimagLite::EnableIf<PsimagLite::IsComplexNumber<ComplexType>::True,ComplexType>::Type
realOrComplexHopping()
{
	return ComplexType(-1,0.1);
}


template<typename SparseMatrixType,typename RealType,typename RngType>
void fillRandomMatrix(SparseMatrixType& t,const RealType& range2, RngType& rng)
{
	typedef typename SparseMatrixType::value_type RealOrComplexType;
	RealType scale = 2.0 + fabs (range2);
	RealOrComplexType hopping = realOrComplexHopping<RealOrComplexType>() / scale;
	SizeType rank = t.rank();

	RealType range = 2.0*range2 / scale;
	SizeType j = 0;
	for (SizeType i = 0; i < rank; i++) {
		t.setRow(i,j);
		t.pushCol(i);
		
		RealType epsilon = range * (rng() - 0.5);
		t.pushValue(epsilon);
		j++;

		SizeType col =  (i + 1) % rank;
		RealOrComplexType hopping2 = (i > col) ? hopping : std::conj(hopping);
		t.pushCol(col);
		t.pushValue(hopping2);
		j++;

		col = (i + rank - 1) % rank;
		hopping2 = (i > col) ? hopping : std::conj(hopping);
		t.pushCol(col);
		t.pushValue(hopping2);
		j++;
	}

	t.setRow(rank,j);
	isHermitian(t,true);
}

template<typename TpemType>
typename TpemType::RealType tpemApply(typename TpemType::TpemSparseType& matrix,
                                      const typename TpemType::BaseFunctorType& funcptr,
                                      const TpemType& tpem)
{
	typedef typename TpemType::RealType RealType;
	SizeType cutoff=tpem.tpemParameters().cutoff;
	std::vector<RealType> coeffs(cutoff);
	std::vector<RealType> moment(cutoff);
	
	tpem.calcCoeffs(coeffs, funcptr);
	tpem.calcMoments(matrix, moment);
	return tpem.expand(moment, coeffs);
}


/* Calculates an observable defined by the function funcptr using exact diagonalization
    for the model given by the Hamiltonian matrix *matrix */
template<typename SparseMatrixType,typename FunctorType>
typename FunctorType::RealType diagApply(const SparseMatrixType& matrix,
                                         const FunctorType& functor)
{
	typedef typename SparseMatrixType::value_type RealOrComplexType;
	typedef typename FunctorType::RealType RealType;
	PsimagLite::Matrix<RealOrComplexType> a;
	crsMatrixToFullMatrix(a,matrix);
	
	std::vector<RealType> eigs(a.n_row());
	diag(a,eigs,'N');

	RealType ret = 0.0;
	for (SizeType i = 0; i < eigs.size(); i++) {
		ret += functor(eigs[i]);
	}	
	return ret;
}

template<typename TpemType,typename FunctorType>
typename FunctorType::RealType tpemApplyDiff(
                     const typename TpemType::TpemSparseType& matrix0,
                     const typename TpemType::TpemSparseType& matrix1,
                     const FunctorType& functor,
                     TpemType& tpem)
{
	typedef typename TpemType::RealType RealType;
	SizeType cutoff=tpem.tpemParameters().cutoff;
	std::vector<RealType> coeffs(cutoff);
	std::vector<RealType> moment(cutoff);

	tpem.calcCoeffs(coeffs, functor);

	tpem.calcMomentsDiff(moment,matrix0,matrix1);
	return tpem.expand(moment, coeffs);
}

template<typename IoInType,typename RealType_>
struct MuBetaStruct {
	typedef RealType_ RealType;
	
	MuBetaStruct(IoInType& io)
	{
		io.rewind();
		io.readline(mu,"Mu=");
		io.readline(beta,"Beta=");
		io.rewind();
	}
	
	MuBetaStruct(const RealType& mu1,const RealType& beta1)
	: mu(mu1),beta(beta1)
	{}

	RealType mu,beta;
};

int main (int argc,char *argv[])
{
	typedef double RealType;
#ifdef USE_COMPLEX
	typedef std::complex<RealType> RealOrComplexType;
#else
	typedef RealType RealOrComplexType;
#endif
	typedef PsimagLite::IoSimple::In IoInType;
	typedef MuBetaStruct<IoInType,RealType> MuBetaStructType;
	typedef Tpem::TpemParameters<IoInType,RealType> TpemParametersType;
	typedef Tpem::Tpem<TpemParametersType,RealOrComplexType> TpemType;
	typedef PsimagLite::CrsMatrix<RealOrComplexType> SparseMatrixType;

	if (argc<2) throw std::runtime_error("Needs the filename\n");

	PsimagLite::Concurrency concurrency(&argc,&argv,1);
	
	IoInType io(argv[1]);
	
	PsimagLite::RandomForTests<RealType> rng(1);
	SparseMatrixType matrix0(400,400);
	fillRandomMatrix(matrix0,10.0,rng);
	
	SparseMatrixType matrix1(400,400);
	fillRandomMatrix(matrix1,10.0,rng);

	MuBetaStructType muBeta(io);
	TpemParametersType tpemParameters(io,muBeta.mu,muBeta.beta);
	std::vector<SizeType> support(2,0);
	support[0] = 0;
	support[1] = matrix0.rank() / 2 - 1;
	tpemParameters.support = support;

	TpemType tpem(tpemParameters);
	SizeType maxCutoff = tpemParameters.cutoff;

// 	tpemOptions.epsProd=1e-5;
// 	tpemOptions.epsTrace=1e-7;
// 	tpemOptions.tpemType="tpem";
// 	tpemOptions.coeffs = "accurate";
	
	std::cout<<"********************************************************\n";
	std::cout<<"****** TESTING TRUNCATED POLYNOMIAL EXPANSION **********\n";
	std::cout<<"********************************************************\n";
	std::cout<<"\n";
	std::cout<<"\n";
	std::cout<<"this testing program calculates model properties in two ways:\n";
	std::cout<<"(i)  Using standard diagonalization\n";
	std::cout<<"(ii) Using the truncated polynomial expansion method\n";
	std::cout<<"\n";
	std::cout<<"All tests are done for a nearest neighbour interaction with\n";
	std::cout<<"random (diagonal) potentials.\n";
	std::cout<<"\n";
	std::cout<<"\n";
	
	matrix1 = matrix0;
	matrix1.setValues(support[0], 2.4 * (rng() - 0.5));
	matrix1.setValues(matrix1.getRowPtr(support[1]), 2.4 * (rng() - 0.5));

	std::cout<<"-------------------------------------------------------------\n";
	std::cout<<"TEST 1: MEAN VALUE FOR THE FUNCTION:                         \n";
	std::cout<<"        E(x) = 5.0 * x * (1.0 - tanh (10.0 * x))\n";
	std::cout<<"\n";

	Tpem::EnergyFunctor<TpemParametersType> energyFunctor(tpemParameters);
	RealType naive0 = diagApply(matrix0, energyFunctor);

	std::cout<<"** Using diagonalization <E>="<<naive0<<"\n";
	std::cout<<"** Using TPEM <E>=(cutoff--> infinity) "
			"lim<E_cutoff> where <E_cutoff> is\n";
	std::cout<<"cutoff\t<E_cutoff>\tError(compared to diag.)\n";

	for (SizeType cutoff = 20; cutoff <= maxCutoff; cutoff++) {
		tpemParameters.cutoff=cutoff;
		RealType tpem0 = tpemApply(matrix0, energyFunctor, tpem);
		std::cout<<cutoff<<"\t"<<tpem0<<"\t"<<naive0<<"\t";
		std::cout<<(100.0 * fabs (1.0 - tpem0 / naive0))<<"\n";
	}

	std::cout<<"-------------------------------------------------------------\n";
	std::cout<<"TEST 2: MEAN VALUE FOR THE FUNCTION:                         \n";
	std::cout<<"        N(x) =  0.5 * (1.0 - tanh (10.0 * x))";
	std::cout<<"\n";

	Tpem::NumberFunctor<TpemParametersType> numberFunctor(tpemParameters);
	naive0 = diagApply(matrix0, numberFunctor);

	std::cout<<"** Using diagonalization <N>="<<naive0<<"\n";
	std::cout<<"** Using TPEM <N>=(cutoff--> infinity) "
			"lim<N_cutoff> where <N_cutoff> is \n";
	std::cout<<"cutoff\t<N_cutoff>\tError(compared to diag.)\n";

	for (SizeType cutoff = 20; cutoff <= maxCutoff; cutoff++) {
		tpemParameters.cutoff=cutoff;
		RealType tpem0 = tpemApply(matrix0, numberFunctor, tpem);
		std::cout<<cutoff<<"\t"<<tpem0<<"\t";
		std::cout<<(100.0 * fabs (1.0 - tpem0 / naive0))<<"\n";
	}
	
	std::cout<<"-------------------------------------------------------------\n";
	std::cout<<"TEST 3: MEAN VALUE AND DIFFERENCE FOR THE FUNCTION:          \n";
	std::cout<<"        S(x) =  log (1.0 + exp (-20.0 * x))\n";
	std::cout<<"\n";
	
	Tpem::ActionFunctor<TpemParametersType> actionFunctor(tpemParameters);
	naive0 = diagApply(matrix0, actionFunctor);
	RealType naive1 = diagApply(matrix1, actionFunctor);

	std::cout<<"** Using diagonalization <S[matrix0]>= "<<naive0<<"\n";
	std::cout<<"** Using diagonalization <S[matrix1]>= "<<naive1<<"\n";
	std::cout<<"** Using diagonalization <S[matrix1]>-<S[matrix0]>="<<(naive0 - naive1)<<"\n";
	std::cout<<"** Using TPEM <S>=(cutoff--> infinity) lim<S_cutoff>\n";
	std::cout<<"cutoff\tDelta_S_cutoff\tS_cutoff[diff]\tError (to diag.)\n";

	for (SizeType cutoff = 10; cutoff <= maxCutoff; cutoff++) {
		tpemParameters.cutoff=cutoff;
		RealType tpem0 = tpemApply(matrix0, actionFunctor, tpem);
		RealType tpem1 = tpemApply(matrix1, actionFunctor, tpem);
		RealType tpemd = tpemApplyDiff(matrix0, matrix1, actionFunctor,tpem);

		std::cout<<cutoff<<"\t"<<(tpem0-tpem1)<<"\t"<<tpemd<<"\t";
		std::cout<<(100.0 * fabs (1.0 - tpemd / (naive0-naive1)))<<"\n";
	}
}
