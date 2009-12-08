/*******************************************************************
 **         THIS FILE IS PART OF THE SPF (SPIN-PHONON-FERMION)    ** 
 **                          COMPUTER PROGRAM                     **
 **                                                               **
 **			 SPF VERSION 6.4                          **
 **                    Dated: March 29, 2006                      **
 **                                                               **
 **   
                                                             **
FOR INTERNAL USE ONLY. 
YOU MAY NOT DISTRIBUTE OR RE-DISTRIBUTE THIS SOFTWARE IN
ANY FORM.
      
DISCLAIMER

THIS SOFTWARE IS PROVIDED BY THE CONTRIBUTORS "AS IS" AND ANY 
EXPRESS OR IMPLIED WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE CONTRIBUTORS BE LIABLE FOR ANY 
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE 
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE 
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, 
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Note: Please acknowledge SPF if you publish results obtained
with it. You can use the following acknowledgement line:
"This research used the SPF 
computer code (http://mri-fre.ornl.gov/spf)." 
 **                                                               **
 **       For more info see: http://mri-fre.ornl.gov/spf          **
 ******************************************************************/




/***************************************************************************************************
* TpemPlus Library by G. A.
* 
****************************************************************************************************/

#ifndef TPEM_H
#define TPEM_H

#ifdef USE_MPI
#include <mpi.h>
#endif

#include <iostream>
#include <string>
#include <vector>
#include <unistd.h>
#include <cstdlib>
#include <cstring>


using std::cerr;
using std::endl;
using std::vector;

#include "tpemplusTypes.h"



 /*!\brief  Structure to store options related to the TPEM.
 
 * Descritpion: Structure that contains TPEMPLUS options.
 * - rank: Rank of the MPI Process. If serial it must be zero.
 * - mpi_nop: Total number of MPI PRocesses. If serial it must be one.
 * - cutoff: Expansion cutoff, M.
 * - epsTrace, epsProd: Expansion epsilons for trace and product truncation.
 * - tpemType: "tpem" for TPEM or "pem" for PEM
 * - coeffs: "alt" for alternative fast or "gsl" for accurate coeffs. with gsl. */
struct TpemOptions
{
        int tpem_rank,rank,mpi_nop1,mpi_nop2;
        unsigned int cutoff;
        double epsTrace,epsProd;
        std::string tpemType,coeffs;
#ifdef USE_MPI
        MPI_Comm mpiCommTpem,mpiCommAdj;
#endif
};

typedef	double	(*tpem_func_ptr)	(double);
double tpem_chebyshev(int m,double x);


tpem_sparse	*new_tpem_sparse	(unsigned int, unsigned int);
void		tpem_sparse_free	(tpem_sparse *);
void		tpem_sparse_copy	(tpem_sparse *, tpem_sparse *);
tpem_t		*tpem_sparse_element	(tpem_sparse *, unsigned int, unsigned int);

void		tpem_matrix_free	(tpem_t **);
//tpem_t		**new_tpem_matrix	(unsigned int, unsigned int);


/*! \brief Holds the subspace in the TPEMPLUS library.
 *
 * struct tpem_subspace is opaque to the user, don't worry about it. But if
 * you must know, the boolean *flags vector tells whether the i^th state has
 * been used so far, think of the macro
 *	#define	isused(i)	(flags[i] != 0)
 * *stack and *top implement a stack of all states that have been used since
 * the last call to tpem_subspace_reset. The rest should be self-explanatory.
 *
 */
typedef	struct	{
	unsigned int	size;
	int	*flags;
	unsigned int	*stack, *top;
} tpem_subspace;


tpem_subspace	*new_tpem_subspace	(unsigned int);
void		tpem_subspace_free	(tpem_subspace *);
void		tpem_subspace_reset	(tpem_subspace *);
void		tpem_subspace_push	(tpem_subspace *, unsigned int);
void tpem_subspace_fill(tpem_subspace *t);

/* tpemplus.cpp */
/*! \brief Initializes the TPEMPLUS library.

* \param argc, argv[]: from int main (INPUT)
* \param tpemOptions: see struct TpemOptions above. (INPUT/OUTPUT)
* \param sfile: ensures that all processes copy argv[1] to sfile. It is optional. (INPUT) */
void tpem_init(int argc,char *argv[],TpemOptions &tpemOptions,char *sfile=0);

/*! \brief Calculates <i|Tm(H)|ket> for all i and all m and store it in moment[i][m]. (EXPERIMENTAL)

* \param matrix: Matrix in sparse format. See struct tpem_sparse above.(INPUT)
* \param moment: A vector of vectors to contain <i|T_m(H)|ket> as moment[i][m].(OUTPUT)
* \param ket: See description.(INPUT)
* \param tpemOptions: see struct TpemOptions above.(INPUT) */
void tpem_off_diagonal_element_tpem (tpem_sparse *matrix, vector<vector<tpem_t> > &moment,
                unsigned int ket, TpemOptions const &tpemOptions);

/*! \brief  Calculates the moments for the expansion.
 * \param matrix: Matrix in sparse format. See struct tpem_sparse above.(INPUT)
 * \param moment: vector<double> that will be filled with the moments.(OUTPUT)
 * \param tpemOptions: see struct TpemOptions above.(INPUT) */
void tpem_calculate_moment (tpem_sparse *matrix, vector<double>  &moment,TpemOptions const &tpemOptions);

/*! \brief Calculates the difference in moments between two matrices that differ only in certain places.
 * \param matrix0, matrix1: The sparse matrices.  See struct tpem_sparse above. (INPUT)
 * \param moment: The difference in moments that will be filled by this function.(OUTPUT)
 * \param support: A list of indices where the matrices matrix0 and matrix1 differ. (INPUT)
 * \param tpemOptions: see struct TpemOptions above.(INPUT) */
void tpem_calculate_moment_diff (tpem_sparse *matrix0, tpem_sparse *matrix1,
		 vector<double> &moment,  vector<unsigned int> const &support,TpemOptions const &tpemOptions);

/* \brief Calculates "fixed" coefficients for the expansion.

 * \param c: coefficients. (OUTPUT)
 * \param f: pointer to function that takes a double and resturns a double.(INPUT) 
 * \param tpemOptions: see struct TpemOptions above.(INPUT) */
void	tpem_calculate_coeffs	(vector<double> &c,tpem_func_ptr f,TpemOptions const &tpemOptions);

/*! \brief Multiplies the moments times the coefficients.

 * \param moments: moments obtained with tpem_calculate_moment (INPUT)
 * \param coeffs: coefficients obtained with tpem_calculate_coeffs (INPUT) */
double tpem_expansion (vector<double> const & moments, vector<double> const &coeffs);

/*! \brief Same as tpem_expansion but for complex moments. */
double tpem_expansion_complex (vector<tpem_t> const &moments, vector<double> const &coeffs);

/*! \brief  It should be called to finalize the TPEMPLUS library. The call is required under MPI. */
void    tpem_finalize();

/*! \brief Broadcast

* On an MPI environment broadcasts dS to all threads and sincronizes all threads by calling MPI_Barrier()
* On a serial environment does nothing.*/
void tpem_bcast(double *dS,TpemOptions const &tpemOptions);
void tpem_bcast(int *dS,TpemOptions const &tpemOptions);
void tpem_bcast(vector<double> &v,TpemOptions const &tpemOptions);

/*
 * 	CALL TREE:
 * 	----------
 *
 *	tpem_calculate_moment:
 *		\tpem_diagonal_element:
 *			\tpem_sparse_product
 *
 *	tpem_calculate_moment_diff:
 *		|tpem_subspace_for_trace:
 *		|	\tpem_diagonal_element:
 *		|		\tpem_sparse_product
 *		\tpem_diagonal_element:
 *			\tpem_sparse_product
 *
 *	tpem_finalize:
 *		|tpem_diagonal_element
 *		|tpem_calculate_moment
 *		\tpem_calculate_moment_diff
 */




#endif
