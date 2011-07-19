
/** \ingroup SPF */
/*@{*/

/*! \file GreenFunction.h
 *
 *  FIXME
 *
 */
#ifndef GREENFUNCTION_H
#define GREENFUNCTION_H
#include "Matrix.h" // in PsimagLite
#include "Fermi.h" // in PsimagLite

namespace Spf {
	template<typename EngineParametersType,typename AlgorithmType>
	class GreenFunction {

	public:
		typedef typename AlgorithmType::FieldType FieldType;
		typedef typename AlgorithmType::ComplexType ComplexType;

		GreenFunction(const EngineParametersType& engineParams,AlgorithmType& algorithm,size_t hilbertSize) :
			engineParams_(engineParams),algorithm_(algorithm),hilbertSize_(hilbertSize),data_(hilbertSize,hilbertSize)
		{
			algorithm_.prepare();
			size_t n = hilbertSize_;
			for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++)
				data_(i,j) = greenFunction(i,j);
		}
		
		ComplexType operator()(size_t lambda1,size_t lambda2)
		{
			return data_(lambda1,lambda2);
		}
		
		FieldType calcNumber() const
		{
			FieldType sum=0;
			for (size_t i=0;i<hilbertSize_;i++) {
				sum += PsimagLite::fermi(
					(algorithm_.e(i)-engineParams_.mu)*engineParams_.beta);
			}
			return sum;
		}
		
		FieldType calcElectronicEnergy() const
		{
			FieldType sum=0;
			for (size_t i=0;i<hilbertSize_;i++) {
				sum +=algorithm_.e(i)*PsimagLite::fermi(
					(algorithm_.e(i)-engineParams_.mu)*engineParams_.beta);
			}
			return sum;
				
		}
		
		void localCharge(std::vector<FieldType>& lc)
		{
			//checkUs();
			//checkLevels();
			for (size_t i=0;i<hilbertSize_;i++) {
				for (size_t lambda=0;lambda<hilbertSize_;lambda++) {
					ComplexType tmp =conj(algorithm_.matrix(i,lambda))*algorithm_.matrix(i,lambda);
					//if (algorithm_.e(lambda)>=engineParams_.mu) continue; // temperature zero
					FieldType s = real(tmp)*PsimagLite::fermi(engineParams_.beta*
							(algorithm_.e(lambda)-engineParams_.mu));
					//if (ether.isSet("savelcd")) {
					//	lc[i+alpha*linSize] = s;
					//} else {
						lc[i] += s;
					//}
				}
				//std::cerr<<"Done i="<<i<<" ss="<<ss<<"\n";
			}
		}

		void electronSpin(std::vector<ComplexType>& es,
				size_t orbitals,size_t n) const
		{
			enum {SPIN_UP,SPIN_DOWN};

			for (size_t i=0;i<n;i++) {
				for (size_t orb=0;orb<orbitals;orb++) {
					size_t x = i+(orb+SPIN_UP*orbitals)*n;
					size_t y = i+(orb+SPIN_DOWN*orbitals)*n;
					// x-component
					es[0] -= data_(x,y) + data_(y,x);
					// y-component
					es[1] += data_(x,y) - data_(y,x);
					// z-component
					es[2] += (1.0 - data_(x,x));
					es[2] -= (1.0 - data_(y,y));
				}
			}
			// y component needs multiplication by sqrt(-1):
			es[1] *= ComplexType(0,1);
		}

		const ComplexType& matrix(size_t lambda1,size_t lambda2) const
		{
			return algorithm_.matrix(lambda1,lambda2);
		}

		const FieldType& e(size_t i) const
		{
			return algorithm_.e(i);
		}

		size_t hilbertSize() const { return algorithm_.hilbertSize(); }
		
		void printMatrix(size_t mode) const
		{
			algorithm_.printMatrix(mode);
		}

	private:
		
		ComplexType greenFunction(size_t lambda1,size_t lambda2) const
		{
			ComplexType sum = 0;
			FieldType beta = engineParams_.beta;
			FieldType mu = engineParams_.mu;

			for (size_t lambda=0;lambda<hilbertSize_;lambda++)
				sum += std::conj(algorithm_.matrix(lambda1,lambda)) *
					algorithm_.matrix(lambda2,lambda) *
				PsimagLite::fermi(-beta*(algorithm_.e(lambda)-mu));
			return sum;
		}

		void checkUs() const
		{
			for (size_t i=0;i<hilbertSize_;i++) {
				ComplexType s = 0;
				for (size_t lambda=0;lambda<hilbertSize_;lambda++) {
					ComplexType tmp =conj(algorithm_.matrix(lambda,i))*algorithm_.matrix(lambda,i);
					s += tmp;
				}
				std::cerr<<"i="<<i<<" sum_lambda U*_{lambda,i} U_{lambda,i} = "<<s<<"\n";
			}
		}

		void checkLevels() const
		{
			FieldType sum = 0;
			for (size_t lambda=0;lambda<hilbertSize_;lambda++) {
				if (algorithm_.e(lambda)>=engineParams_.mu) continue;
				sum++;
			}
			std::cerr<<"mu = "<<engineParams_.mu<<" Levels below="<<sum<<"\n"; 
		}

		const EngineParametersType& engineParams_;	
		AlgorithmType& algorithm_;
		size_t hilbertSize_;
		PsimagLite::Matrix<ComplexType> data_;
		
	}; // GreenFunction
	
	template<typename EngineParametersType,typename AlgorithmType>
	std::ostream& operator<<(std::ostream& os,GreenFunction<EngineParametersType,AlgorithmType>& algorithm)
	{
		throw std::runtime_error("unimplemented operator<< for GreenFunction\n");
		return os;
	}
} // namespace Spf

/*@}*/
#endif
