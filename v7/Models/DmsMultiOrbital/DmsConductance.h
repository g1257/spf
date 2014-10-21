/** \ingroup SPF */
/*@{*/

/*! \file DmsConductance.h
 *
 *
 *
 */
#ifndef DMS_CONDUCTANCE_H
#define DMS_CONDUCTANCE_H
#include "Vector.h"

namespace Spf {

template<typename MatrixType, typename GreenFunctionType>
class DmsConductance {

	typedef typename MatrixType::value_type ComplexType;
	typedef typename ComplexType::value_type RealType;
	typedef PsimagLite::Vector<int>::Type VectorIntType;

public:

	RealType operator()(const GreenFunctionType& gf,
	                    RealType mu,
	                    int dim,
	                    int l1,
	                    int maxiter,
	                    RealType maxerror)
	{
		int n = gf.hilbertSize();
		VectorIntType pivot(n);
		MatrixType h0(n,n);
		MatrixType h1(n,n);
		MatrixType t0(n,n);
		MatrixType t1(n,n);
		MatrixType x(n,n);
		MatrixType sigmaRight(n,n);
		MatrixType sigmaLeft(n,n);
		MatrixType greenPlus(n,n);

		/* intralayer (resp. interlayer) hamiltonian h0 (resp. h1) */

		switch (dim) {
		case 3:
			cond3d(x,h0,h1,gf,l1);
			break;
		case 2:
			cond2d(x,h0,h1,gf,l1);
			break;
		default:
			throw PsimagLite::RuntimeError("Cannot calculate conductance in 1d\n");
		}

		/* velocity operator */
		ComplexType one = 1.0;
		ComplexType zero = 0.0;
		char nChar = 'n';
		psimag::BLAS::zgemm_(&nChar, &nChar, &n, &n, &n,&one, &(h0(0,0)), &n,
		                     &(x(0,0)), &n,&zero, &(t0(0,0)), &n);
		psimag::BLAS::zgemm_(&nChar, &nChar, &n, &n, &n,&one, &x(0,0), &n,
		                     &h0(0,0), &n,&zero, &t1(0,0), &n);

		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; j++) {
				x(i,j) = (t0(i,j) - t1(i,j)) /(l1 - 1.0);
			}
		}

		RealType eta = 0.0;
		ComplexType e(mu, eta);

		/* right self-energy */
		for (int i = 0; i < n; ++i)
			sigmaRight(i,i) = ComplexType(0.0, -1.0);

		RealType error = 1.0;
		int m = n * n;
		char cChar = 'c';
		for (int iter = 0; iter < maxiter && error > maxerror; ++iter) {

			for (int i = 0; i < n; ++i) {
				for (int j = 0; j < n; ++j) {
					if (i!=j)
						t1(i,j) =  0.0- h0(i,j)- sigmaRight(i,j);
					else
						t1(i,j) = e - h0(i,j)- sigmaRight(i,j);

				}
			}

			int x0 = 0;
			psimag::LAPACK::zgetrf_(&n, &n, &t1(0,0), &n, &pivot[0], &x0);
			psimag::LAPACK::zgetri_(&n, &t1(0,0), &n, &pivot[0], &t0(0,0), &m, &x0);
			/* t1 = h1 * t1 * h1^H */
			psimag::BLAS::zgemm_(&nChar, &nChar, &n, &n, &n, &one, &t1(0,0), &n,
			                     &h1(0,0), &n,&zero, &t0(0,0), &n);
			psimag::BLAS::zgemm_(&cChar, &nChar, &n, &n, &n, &one, &h1(0,0), &n,
			                     &t0(0,0), &n,&zero, &t1(0,0), &n);

			error = 0.0;
			for (int i = 0; i < n; ++i) {
				for (int j = 0; j < n; ++j) {
					t1(i,j) = 0.5*(t1(i,j) + sigmaRight(i,j));

					error += std::norm(t1(i,j) - sigmaRight(i,j));
					sigmaRight(i,j) = t1(i,j);
				}
			}

			error /= n * n;
		}

		/* left self-energy */
		for (int i = 0; i < n; ++i)
			sigmaLeft(i,i) = ComplexType(0.0, -1.0);

		error = 1.0;
		for (int iter = 0; iter < maxiter && error > maxerror; ++iter) {

			for (int i = 0; i < n; ++i) {
				for (int j = 0; j < n; ++j) {
					if (i!=j)
						t1(i,j) = 0.0 - h0(i,j);
					else
						t1(i,j) = e - h0(i,j);
					t1(i,j) = t1(i,j)- sigmaLeft(i,j);

				}
			}

			int x0 = 0;
			psimag::LAPACK::zgetrf_(&n, &n, &t1(0,0), &n, &pivot[0], &x0);
			psimag::LAPACK::zgetri_(&n, &t1(0,0), &n, &pivot[0], &t0(0,0), &m, &x0);
			psimag::BLAS::zgemm_(&nChar, &cChar, &n, &n, &n, &one, &t1(0,0), &n,
			                     &h1(0,0), &n,&zero, &t0(0,0), &n);
			psimag::BLAS::zgemm_(&nChar, &nChar, &n, &n, &n, &one, &h1(0,0), &n,
			                     &t0(0,0), &n, &zero, &t1(0,0), &n);

			error = 0.0;
			for (int i = 0; i < n; ++i) {
				for (int j = 0; j < n; ++j) {
					t1(i,j) = 0.5*(t1(i,j) + sigmaLeft(i,j));

					error += std::norm(t1(i,j) - sigmaLeft(i,j));
					sigmaLeft(i,j) = t1(i,j);
				}
			}

			error /= n * n;
		}

		/* construct G+ */
		for (int i = 0; i < n; ++i)
			greenPlus(i,i) = ComplexType(mu, 0.0);

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				greenPlus(i,j)=greenPlus(i,j)-h0(i,j);
			}
		}

		/* self energies */
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				greenPlus(i,j)=greenPlus(i,j)-sigmaRight(i,j);
				greenPlus(i,j)=greenPlus(i,j)-sigmaLeft(i,j);

			}
		}

		int x0 = 0;
		psimag::LAPACK::zgetrf_(&n, &n, &greenPlus(0,0), &n, &pivot[0], &x0);
		if (x0) {
			PsimagLite::String str("FATAL -- INTERNAL --\n"
			                       "Problems factorizing greenPlus\n");
			throw PsimagLite::RuntimeError(str);
		}

		psimag::LAPACK::zgetri_(&n, &greenPlus(0,0), &n, &pivot[0], &t0(0,0), &m, &x0);
		if (x0) {
			PsimagLite::String str("FATAL -- INTERNAL --\n"
			                       "Problems inverting greenPlus\n");
		}

		/* construct G+ - G- */
		for (int i = 0; i < n; i++) {
			for (int j = 0; j <= i; j++) {
				error= 0.5 * (imag(greenPlus(j,i))
				              + imag(greenPlus(i,j)));
				greenPlus(j,i) = ComplexType(error,0.5 * (-real(greenPlus(j,i))
				                                          + real(greenPlus(i,j))));

			}
		}

		for (int i = 0; i < n; i++) {
			for (int j = 0; j <= i; j++) {
				greenPlus(i,j) =  ComplexType(real(greenPlus(j,i)),
				                              -imag(greenPlus(j,i)));
			}
		}

		/* Kubo formula */
		psimag::BLAS::zgemm_(&nChar, &nChar, &n, &n, &n, &one, &x(0,0), &n,
		                     &greenPlus(0,0), &n,&zero, &t0(0,0),  &n);
		psimag::BLAS::zgemm_(&nChar, &nChar, &n, &n, &n, &one, &t0(0,0), &n,
		                     &t0(0,0), &n, &zero, &greenPlus(0,0), &n);

		RealType cond = 0.0;
		for (int i = 0; i < n; i++) {
			cond +=  real(greenPlus(i,i));
		}

		return (2.0 *cond);
	}

private:

	void cond3d(MatrixType& x,
	            MatrixType& h0,
	            MatrixType& h1,
	            const GreenFunctionType& matrix,
	            SizeType l1)
	{
		SizeType n = x.n_row();
		SizeType volume = n/2;

		for (SizeType i = 0; i < l1; i++) {
			for (SizeType j = 0; j < l1; j++) {
				for (SizeType k=0;k<l1;k++) {
					SizeType x0=k*l1*l1+j*l1+i;
					x(x0,x0)=x(x0+volume,x0+volume)=k;
				}
			}
		}

		for (SizeType i = 0; i < l1; i++) {
			for (SizeType j = 0; j < l1; j++) {
				for (SizeType k = 0; k < l1; k++) {
					SizeType x1 = k*l1*l1 + j*l1 +i;
					for (SizeType x0=0;x0<l1;x0++) {
						for (SizeType y0=0;y0<l1;y0++) {
							for (SizeType z0=0;z0<l1;z0++) {
								SizeType y1=z0*l1*l1+y0*l1+x0;
								if (z0==0 && k==l1-1) {
									h1(x1,y1)=matrix(x1,y1);
									h1(x1+volume,y1+volume)=matrix(x1+volume,y1+volume);
									h1(x1+volume,y1)=matrix(x1+volume,y1);
									h1(x1,y1+volume)=matrix(x1,y1+volume);
								} else if (k!=0 || z0!=l1-1) {
									h0(x1,y1)=matrix(x1,y1);
									h0(x1+volume,y1+volume)=matrix(x1+volume,y1+volume);
									h0(x1+volume,y1)=matrix(x1+volume,y1);
									h0(x1,y1+volume)=matrix(x1,y1+volume);

								}
							}
						}
					}
				}
			}
		}
	}

	void cond2d(MatrixType& x,
	            MatrixType& h0,
	            MatrixType& h1,
	            const GreenFunctionType& matrix,
	            SizeType l1)
	{
		SizeType n = x.n_row();
		SizeType volume = n/2;

		for (SizeType i = 0; i < l1; i++) {
			for (SizeType j = 0; j < l1; j++) {
				SizeType x0=j*l1+i;
				x(x0,x0)=x(x0+volume,x0+volume)=j;
			}
		}

		for (SizeType i = 0; i < l1; i++) {
			for (SizeType j = 0; j < l1; j++) {
				SizeType x1 =j*l1 +i;
				for (SizeType x0=0;x0<l1;x0++) {
					for (SizeType y0=0;y0<l1;y0++) {
						SizeType y1=y0*l1+x0;
						if (y0==0 && j==l1-1) {
							h1(x1,y1)=matrix(x1,y1);
							h1(x1+volume,y1+volume)=matrix(x1+volume,y1+volume);
							h1(x1+volume,y1)=matrix(x1+volume,y1);
							h1(x1,y1+volume)=matrix(x1,y1+volume);
						} else if (j!=0 || y0!=l1-1) {
							h0(x1,y1)=matrix(x1,y1);
							h0(x1+volume,y1+volume)=matrix(x1+volume,y1+volume);
							h0(x1+volume,y1)=matrix(x1+volume,y1);
							h0(x1,y1+volume)=matrix(x1,y1+volume);

						}
					}
				}
			}
		}
	}
}; // class DmsConductance

} // namespace Spf

#endif // DMS_CONDUCTANCE_H

