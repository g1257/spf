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



#include "tpemplus.h"

extern "C" {
#ifndef NO_GSL
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#endif
}

struct my_f_params { int m;  tpem_func_ptr funk;};


void my_handler (const char * reason, const char * file, int line, int gsl_errno);

#ifdef USE_MPI
MPI_Comm adj;

void MPI_VecRecv(vector<double> &v, int n, MPI_Datatype type1,int rank, int tag,MPI_Comm comm,MPI_Status *status )
{
	double *vv = new double[n];
	
	MPI_Recv(vv,n,type1,rank,tag,comm,status);
	for (unsigned int i=0;i<n;i++) v[i]=vv[i];
	
	delete [] vv;
}

void MPI_VecSend(vector<double> &v, int n, MPI_Datatype type1, int rank,int tag,MPI_Comm comm)
{
	double *vv = new double[n];
	for (unsigned int i=0;i<n;i++) vv[i]=v[i];
	MPI_Send(vv,n,type1,rank,tag,comm);
	delete [] vv;
}

void MPI_VecReduce(vector<double> const &v1,vector<double> &v2,int n,MPI_Datatype type1,MPI_Op op1,int rank,MPI_Comm comm)
{

	double *vv1 = new double[n];
	double *vv2 = new double[n];
	for (unsigned int i=0;i<n;i++) vv1[i]=v1[i];
	MPI_Reduce(vv1,vv2,n,type1,op1,rank,comm);
	for (unsigned int i=0;i<n;i++) v2[i]=vv2[i];
	delete [] vv1;
	delete [] vv2;
}

void MPI_VecReduce(vector<tpem_t> const &v1,vector<tpem_t> &v2,int n,MPI_Datatype type1,MPI_Op op1,int rank,MPI_Comm comm)
{

	double *vv1r = new double[n];
	double *vv1i = new double[n];
	double *vv2r = new double[n];
	double *vv2i = new double[n];
	unsigned int i;
	for (i=0;i<n;i++) {
		vv1r[i]=real(v1[i]);
		vv1i[i]=imag(v1[i]);
	}
	MPI_Reduce(vv1r,vv2r,n,type1,op1,rank,comm);
	MPI_Reduce(vv1i,vv2i,n,type1,op1,rank,comm);
	for (i=0;i<n;i++) v2[i]=tpem_t(vv2r[i],vv2i[i]);
		
	
	delete [] vv1r;
	delete [] vv2r;
	delete [] vv1i;
	delete [] vv2i;
}

void MPI_VecRecv(vector<tpem_t> &v, int n, char *type1,int rank, int tag,MPI_Comm comm,MPI_Status *status )
{
	unsigned int nn=v.size();
	double *vreal = new double[nn];
	double *vimag = new double[nn];
	int tag1=12,tag2=15;
	
	MPI_Recv(vreal,n,MPI_DOUBLE,rank,tag1,comm,status);
	MPI_Recv(vimag,n,MPI_DOUBLE,rank,tag2,comm,status);
	
	for (unsigned int i=0;i<nn;i++) v[i]=tpem_t(vreal[i],vimag[i]);
	delete [] vreal;
	delete [] vimag;
}

void MPI_VecSend(vector<tpem_t> &v, int n, char *type1, int rank,int tag,MPI_Comm comm)
{
	unsigned int nn=v.size();
	double *vreal = new double[nn];
	double *vimag = new double[nn];
	int tag1=12,tag2=15;
	
	for (unsigned int i=0;i<nn;i++) {
		vreal[i]=real(v[i]);
		vimag[i]=imag(v[i]);
	}
	
	MPI_Send(vreal,n,MPI_DOUBLE,rank,tag1,comm);
	MPI_Send(vimag,n,MPI_DOUBLE,rank,tag2,comm);
	
	delete [] vreal;
	delete [] vimag;
}
void MPI_VecReduce(vector<double> &v1)
{
	MPI_VecReduce(v1,v1,v1.size(),MPI_DOUBLE,MPI_SUM,0,adj);
}

#endif

void tpem_init(int argc, char *argv[], TpemOptions &tpemOptions,char  *sfile) 
{
	//int i;
	int MAX_WORD=256;
	char *tempc;
	int mpi_nop1=1,tpem_rank=0,rank=0;
	//int mpi_nop=1,tmp;
	
	tempc = new char[MAX_WORD]; 
	
	tpemOptions.mpi_nop1=mpi_nop1;
#ifdef USE_MPI
	MPI_Comm newcomm;
	
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&mpi_nop);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	
	if (rank==0) {
		std::cout<<"**************************************************************************\n";
		std::cout<<"This program uses TpemPlus Library v1.0 by G.A.\n";
		std::cout<<"**************************************************************************\n\n";
		if (argc>=3) {
			std::cout<<"tpem_init: rank==0 got argc==3 and so setting up mpi_nop2\n";
			tpemOptions.mpi_nop2 = atoi(argv[2]);
		} else {
			tpemOptions.mpi_nop2 = 1;
		}
	}
	MPI_Bcast(&(tpemOptions.mpi_nop2),1,MPI_INT,0,MPI_COMM_WORLD);
	
	if (mpi_nop % tpemOptions.mpi_nop2 !=0) {
		cerr<<"tpem_init: mpi_nop2="<<tpemOptions.mpi_nop2<<" must divide mpi_nop="<<mpi_nop<<endl;
		MPI_Finalize();
		exit(1);
	}	
	
	tpemOptions.mpi_nop1 = mpi_nop/tpemOptions.mpi_nop2;
	getcwd(tempc,MAX_WORD);
        MPI_Bcast(tempc,strlen(tempc)+1,MPI_CHAR,0,MPI_COMM_WORLD);
        if (rank!=0) chdir(tempc);
	MPI_Barrier(MPI_COMM_WORLD);
	getcwd(tempc,MAX_WORD);	
	MPI_Bcast(&argc,1,MPI_INT,0,MPI_COMM_WORLD);
	
	for (i=0;i<argc;i++) {
		MPI_Bcast(argv[i],strlen(argv[i])+1,MPI_CHAR,0,MPI_COMM_WORLD);
	}
	delete [] tempc;
	
	// setup tpem group
	if (tpemOptions.mpi_nop2>1) {
		int color = int(rank/tpemOptions.mpi_nop1);
		MPI_Comm_split (MPI_COMM_WORLD, color, rank, &newcomm);
		MPI_Comm_rank(newcomm,&tpem_rank);
		MPI_Comm_size(newcomm,&tmp);
		if (tmp!=tpemOptions.mpi_nop1) {
			cerr<<"tpem_init: Internal error: tmp="<<tmp<<"!=mpi_nop2="<<tpemOptions.mpi_nop2<<endl;
			tpem_finalize();
			exit(1);
		}
		
		cerr<<"tpem_init: proc has rank="<<rank<<" tpem_rank="<<tpem_rank<<endl;
		tpemOptions.mpiCommTpem=newcomm;
	
		// setup adjoint group
		if (rank==0 || (rank%tpemOptions.mpi_nop1)==0) {
			color=0;
		} else {
			color=1;
		}
		MPI_Comm_split (MPI_COMM_WORLD, color, rank, &newcomm); 
		adj=tpemOptions.mpiCommAdj=newcomm;
	} else {
		tpemOptions.mpiCommTpem=MPI_COMM_WORLD;
		adj=tpemOptions.mpiCommAdj=MPI_COMM_WORLD;
		tpem_rank = rank;
	}
#else
	if (tpemOptions.mpi_nop2!=1) {
		cerr<<"tpem_init: mpi_nop2="<<tpemOptions.mpi_nop2<<" but 1 expected. Either set mpi_nop2=1 or use MPI\n";
		exit(1);
	}
#endif	
	
	tpemOptions.rank=rank;
	tpemOptions.tpem_rank=tpem_rank;
	if (sfile) strcpy(sfile,argv[1]);
	//std::cerr<<"tpem_init: process "<<rank<<" has read file "<<sfile<<endl;
	
#ifndef NO_GSL
	gsl_set_error_handler (&my_handler);
#endif  	
	
}


void vectorAcc(vector<double> &dest,vector<double> const &src)
{
        int j,cutoff=src.size();
        for (j=0;j<cutoff;j++) {
                dest[j] += src[j];
        }
}

void vectorAcc(vector<tpem_t> &dest,vector<tpem_t> const &src)
{
        int j,cutoff=src.size();
        for (j=0;j<cutoff;j++) {
                dest[j] = dest[j] + src[j];
        }
}



inline double tpem_norm (tpem_t z)
{
#ifdef	TPEM_COMPLEX
	return	real(z)*real(z) + imag(z)*imag(z);
#else
	return	z * z;
#endif
}

inline tpem_t tpem_conj_mult (tpem_t x, tpem_t y)
{
	return conj (x)*y;
}


void tpem_sparse_product_tpem (tpem_sparse *matrix, tpem_subspace *info,
		vector<tpem_t> &dest, vector<tpem_t> const &src, TpemOptions const &tpemOptions,
		double eps)
{
	unsigned int	*oldtop, *p, i, j, k;
	tpem_t	t;
	double	u;

	for (i = 0; i < matrix->rank; i++)
		dest[i] = 0.0;

	oldtop = info->top;

	/* loop over states that have been used so far */
	for (p = info->stack; p < oldtop; p++) {
		j = *p;
		//if (real(src[j])==0 && imag(src[j])==0) continue;
		/* loop over nonzero elements of j^th row */
		for (k = matrix->rowptr[j]; k < matrix->rowptr[j + 1]; k++) {
			i = matrix->colind[k];
			t = src[j] * matrix->values[k];
			u = tpem_norm (t);
			dest[i] += t;
			if (u > eps) tpem_subspace_push (info, i);
		}
	}
}

void tpem_sparse_product_pem (tpem_sparse *matrix,vector<tpem_t> &dest, vector<tpem_t> const &src,
	TpemOptions const &tpemOptions)
{	
	unsigned int	i, j, k;
	tpem_t	t;
	
	for (i = 0; i < matrix->rank; i++)
		dest[i] = 0.0;

	
	/* loop over all rows */
	for (j=0;j<matrix->rank;j++) {
		/* loop over nonzero elements of j^th row */	
		for (k = matrix->rowptr[j]; k < matrix->rowptr[j + 1]; k++) {
#ifdef FAST_MULTIPLIER
                        if (tpem_norm(src[j])<1e-6) continue;
#endif
			i = matrix->colind[k];
			t = src[j] * matrix->values[k];
			dest[i] = dest[i] + t;
		}
	}
}


void tpem_diagonal_element_tpem (tpem_sparse *matrix,
		vector<double> &moment, unsigned int ket, TpemOptions const &tpemOptions,
		double eps,tpem_subspace *info)
{	
	unsigned int n=moment.size();
	vector<tpem_t>	tmp(matrix->rank), jm0(matrix->rank), jm1(matrix->rank);
	unsigned int	i, m, *p;
	tpem_t	sum1, sum2, keep;
	

	for (i = 0; i < matrix->rank; i++)
		tmp[i] = jm0[i] = jm1[i] = 0.0;
	jm0[ket] = 1.0;	/* set |j,0> */
	
	tpem_subspace_reset (info);
	tpem_subspace_push (info, ket);

	/* calculate |j,1> = X|j,0> */
	tpem_sparse_product_tpem (matrix, info, jm1, jm0, tpemOptions,tpemOptions.epsProd);
	
	sum1 = tpem_conj_mult (jm0[ket], jm1[ket]);
	moment[1] += real (sum1);
	
	sum2 = 0.0;
	for (p = info->stack; p < info->top; p++)
		sum2 = sum2 + tpem_conj_mult (jm1[*p], jm1[*p]);
	moment[2] += real (sum2);

	/* calculate |j,m> = 2X|j,m-1> - |j,m-2>
	 *
	 * begin (m=2) pass	jm0 = |j,0>	jm1 = |j,1>
	 * end   (m=2) pass	jm0 = |j,1>	jm1 = |j,2>
	 * ...
	 * begin (m=k) pass	jm0 = |j,k-2>	jm1 = |j,k-1>
	 * end   (m=k) pass	jm0 = |j,k-1>	jm1 = |j,k>
	 */
	for (m = 2; m < n / 2; m++) {
		/* calculate |tmp> = X|jm1> */
		tpem_sparse_product_tpem (matrix, info, tmp, jm1, tpemOptions,eps);
		sum1 = sum2 = 0.0;
		for (p = info->stack; p < info->top; p++) {
			i = *p;
			keep = tmp[i] + tmp[i]  - jm0[i];
			/* for moment[2 * m    ] */
			sum1 = sum1 + tpem_conj_mult (keep, keep);
			/* for moment[2 * m - 1] */
			sum2 = sum2 + tpem_conj_mult (keep, jm1[i]);
			/* set |j,m-1> and |j,m-2> for next iteration */
			jm0[i] = jm1[i];
			jm1[i] = keep;
		}
		moment[m + m    ] += real (sum1);
		moment[m + m - 1] += real (sum2);
		if (fabs(imag(sum2))>1e-6) {
			cerr<<"tpem_diagonal_element_tpem: problem "<<imag(sum2)<<endl;
		}
	}
}

void  tpem_diagonal_element_pem (tpem_sparse *matrix, vector<double> &moment, unsigned int ket,
	TpemOptions const &tpemOptions)
{	
	unsigned int n=moment.size();
	vector<tpem_t>	tmp(matrix->rank), jm0(matrix->rank), jm1(matrix->rank);
	unsigned int	i, j,m;
	tpem_t	sum1, sum2, keep;
	
	for (i = 0; i < matrix->rank; i++)
		tmp[i] = jm0[i] = jm1[i] = 0.0;
	jm0[ket] = 1.0;	/* set |j,0> */
	
	
	/* calculate |j,1> = X|j,0> */
	tpem_sparse_product_pem (matrix, jm1, jm0,tpemOptions);
	
	sum1 = tpem_conj_mult (jm0[ket], jm1[ket]);
	moment[1] += real (sum1);
	
	sum2 =  0.0;
	for (j=0;j<matrix->rank;j++)
		sum2 = sum2 + tpem_conj_mult (jm1[j], jm1[j]);
	moment[2] += real (sum2);

	/* calculate |j,m> = 2X|j,m-1> - |j,m-2>
	 *
	 * begin (m=2) pass	jm0 = |j,0>	jm1 = |j,1>
	 * end   (m=2) pass	jm0 = |j,1>	jm1 = |j,2>
	 * ...
	 * begin (m=k) pass	jm0 = |j,k-2>	jm1 = |j,k-1>
	 * end   (m=k) pass	jm0 = |j,k-1>	jm1 = |j,k>
	 */
	for (m = 2; m < n / 2; m++) {
		/* calculate |tmp> = X|jm1> */
		tpem_sparse_product_pem (matrix, tmp, jm1,tpemOptions);
		sum1 = sum2 = 0.0;
		for (i=0;i<matrix->rank;i++) {
			keep = tmp[i] + tmp[i] - jm0[i];
			/* for moment[2 * m    ] */
			sum1 = sum1 +  tpem_conj_mult (keep, keep);
			/* for moment[2 * m - 1] */
			sum2 = sum2 + tpem_conj_mult (keep, jm1[i]);
			/* set |j,m-1> and |j,m-2> for next iteration */
			jm0[i] = jm1[i];
			jm1[i] = keep;
		}
		moment[m + m    ] += real (sum1);
		moment[m + m - 1] += real (sum2);
	}
}

void tpem_diagonal_element(tpem_sparse *matrix,
		vector<double> &moment, unsigned int ket, TpemOptions const &tpemOptions)
{
	
	
	tpem_subspace *work = new_tpem_subspace (matrix->rank);
	if (tpemOptions.tpemType == "tpem") {
		tpem_diagonal_element_tpem (matrix,moment,ket, tpemOptions,tpemOptions.epsProd,work);
	} else if (tpemOptions.tpemType == "pem") {
		tpem_diagonal_element_pem (matrix,moment,ket, tpemOptions);
	} else {
		cerr<<"tpem_diagonal_element: Unknown type: "<<tpemOptions.tpemType<<endl;
		exit(1);
	}
	tpem_subspace_free(work);
}

void tpem_diagonal_element(tpem_sparse *matrix,
		vector<double> &moment, unsigned int ket, TpemOptions const &tpemOptions,
		tpem_subspace *info)
{
	
	if (!(tpemOptions.tpemType == "tpem")) {
		cerr<<"tpem_diagonal_element: Expected tpemType==tpem but tpemType==";
		cerr<<tpemOptions.tpemType<<endl;
		exit(1);
	}
	tpem_diagonal_element_tpem (matrix,moment,ket, tpemOptions,tpemOptions.epsTrace,info);
}

void tpem_subspace_for_trace (
		tpem_sparse *matrix0, tpem_sparse *matrix1,
		vector<double> &moment0, vector<double> &moment1,
		vector<unsigned int> const  &support, TpemOptions const &tpemOptions,
		tpem_subspace	*info)
{	/* f77: mkTraceState */
	
	
	tpem_subspace	*work;
	unsigned int	i, j, *p;
	unsigned int nsupp=support.size();
	
	work = new_tpem_subspace (matrix0->rank);
	tpem_subspace_reset (info);
	
	for (i = 0; i < nsupp; i++) {
		j = support[i];
		tpem_diagonal_element (matrix0, moment0, j,tpemOptions,work);
		for (p = work->stack; p < work->top; p++)
			tpem_subspace_push (info, *p);
		tpem_diagonal_element (matrix1, moment1, j, tpemOptions,work);
		for (p = work->stack; p < work->top; p++)
			tpem_subspace_push (info, *p);
	}
	
	tpem_subspace_free(work);
}

void tpem_calculate_moment (tpem_sparse *matrix, vector<double>  &moment,TpemOptions const &tpemOptions)
{	
	unsigned int n=moment.size();
	unsigned int	i,k,part,activeProc;
	vector<double> buf(n);
	int rank = tpemOptions.tpem_rank;
	int mpi_nop1 = tpemOptions.mpi_nop1;
	//double dummy=0.0;
	
	
	for (i = 0; i < n; i++)
		moment[i] = buf[i] = 0.0;
	moment[0] = (double) matrix->rank;
	
	part=matrix->rank/mpi_nop1;
	if (matrix->rank % mpi_nop1>0) part++;
	activeProc=matrix->rank/part;
	if (matrix->rank % part>0) activeProc++;
	
	for (k=0;k<part;k++) {
		i=rank*part+k;
		if (i>=matrix->rank) break;
		tpem_diagonal_element (matrix, moment, i, tpemOptions);
	}
	
	
#ifdef USE_MPI
	MPI_Status status;
	int tag=0;
	
	MPI_VecReduce(moment,moment,n,MPI_DOUBLE,MPI_SUM,0,tpemOptions.mpiCommTpem);
	
	moment[0] = (double) matrix->rank;
	MPI_Barrier(tpemOptions.mpiCommTpem);
#endif

	if (rank==0) {for (i = 2; i < n; i += 2)
		moment[i] = 2.0 * moment[i] - moment[0];
	
	for (i = 3; i < n - 1; i += 2)
		moment[i] = 2.0 * moment[i] - moment[1];
	}

}


void tpem_calculate_moment_diff (tpem_sparse *matrix0, tpem_sparse *matrix1,
		 vector<double> &moment,  vector<unsigned int> const &support,TpemOptions const &tpemOptions)
{	
	
	unsigned int n=moment.size();
	vector<double>	moment0(n), moment1(n),buf0(n),buf1(n);
	tpem_subspace	*info;
	unsigned int	i, *p;	
	unsigned int part,k,resto,activeProc;
	int rank=tpemOptions.tpem_rank;
	int mpi_nop1=tpemOptions.mpi_nop1;
	//double dummy=0.0;
	unsigned int total;
	
	info = new_tpem_subspace(matrix0->rank);
	if (tpemOptions.tpemType=="tpem") {
		tpem_subspace_for_trace (matrix0, matrix1, moment0, moment1,support,tpemOptions,info);
	} else {
		tpem_subspace_fill(info);
	}
	
	for (i = 0; i < n; i++)
		moment0[i] = moment1[i] = 0.0;
	moment0[0] = moment1[0] = (double) matrix0->rank;
	
	total = info->top-info->stack;
	if (size_t(mpi_nop1)<2*total) {		
		part  = total/mpi_nop1;
		resto =	total % mpi_nop1;
		if (resto>0) part++;
		activeProc=total/part;
		if (total % part >0) activeProc++;
	} else {
		part=1;
		activeProc=2*total;
	}
	for (k=0;k<part;k++) {
		if (size_t(mpi_nop1)<2*total) {
			p=info->stack + (part*rank) +k;
			if (p>=info->top) {
				moment0[0]=moment1[0]=0.0;
				break;
			}
			tpem_diagonal_element (matrix0, moment0, *p, tpemOptions);
			tpem_diagonal_element (matrix1, moment1, *p, tpemOptions);
		} else {
			if (size_t(rank)<total) p=info->stack + (part*rank) +k;	
			else p=info->stack + (part*(rank-total)) +k;
			if (p>=info->top) {
				moment0[0]=moment1[0]=0.0;
				break;
			}
			if (size_t(rank)<total) tpem_diagonal_element (matrix0, moment0, *p, tpemOptions);
			else tpem_diagonal_element (matrix1, moment1, *p, tpemOptions);
		}
	}
	
#ifdef USE_MPI
	MPI_Status mstatus;
	int tag0=98,tag1=99;
	
	if (mpi_nop1<2*total) {
		MPI_VecReduce(moment0,moment0,n,MPI_DOUBLE,MPI_SUM,0,tpemOptions.mpiCommTpem);
		MPI_VecReduce(moment1,moment1,n,MPI_DOUBLE,MPI_SUM,0,tpemOptions.mpiCommTpem);
	} else {
		// moment0=0 for rank>=total execpt moment0[0]
		MPI_VecReduce(moment0,moment0,n,MPI_DOUBLE,MPI_SUM,0,tpemOptions.mpiCommTpem);
		// moment1=0 for rank<total moment1[0]
		MPI_VecReduce(moment1,moment1,n,MPI_DOUBLE,MPI_SUM,total,tpemOptions.mpiCommTpem);
		if (rank==total) MPI_VecSend(moment1,n,MPI_DOUBLE,0,tag0,tpemOptions.mpiCommTpem);
		if (rank==0) MPI_VecRecv(moment1,n, MPI_DOUBLE,total,tag0,tpemOptions.mpiCommTpem,&mstatus);
	}
	
	MPI_Barrier(tpemOptions.mpiCommTpem);
	moment0[0] = moment1[0] = (double) matrix0->rank;
#endif
	

	for (i = 2; i < n; i += 2) {
		moment0[i] = 2.0 * moment0[i] - moment0[0];
		moment1[i] = 2.0 * moment1[i] - moment1[0];
	}

	for (i = 3; i < n - 1; i += 2) {
		moment0[i] = 2.0 * moment0[i] - moment0[1];
		moment1[i] = 2.0 * moment1[i] - moment1[1];
	}

	for (i = 0; i < n; i++)
		moment[i] = moment0[i] - moment1[i];
		
	tpem_subspace_free(info);	
	
}


double tpem_expansion (vector<double> const & moments, vector<double> const &coeffs) 
{
	unsigned int	i,n=moments.size();
	double	ret = 0.0;
	
	if (n!=coeffs.size()) {
		cerr<<"tpem_expansion: coeffs and moments of different size\n";
	}
	for (i = 0; i < n; i++)
		ret += moments[i] * coeffs[i];
	return ret;
}

double tpem_expansion_complex (vector<tpem_t> const &moments, vector<double> const &coeffs) 
{
	unsigned int	i,n=moments.size();
	double	ret = 0.0;
	
	if (n!=coeffs.size()) {
		cerr<<"tpem_expansion_complex: coeffs and moments of different size\n";
	}
	for (i = 0; i < n; i++)
		ret += real(moments[i]) * coeffs[i];
	return ret;
}
	
void tpem_bcast(vector<double> &v,TpemOptions const &tpemOptions)
{
#ifdef USE_MPI
	unsigned int i;
	double *vv = new double[v.size()];
	for (i=0;i<v.size();i++) vv[i]=v[i];
	MPI_Bcast(vv,v.size(),MPI_DOUBLE,0,tpemOptions.mpiCommTpem);       
	MPI_Barrier(tpemOptions.mpiCommTpem);
	for (i=0;i<v.size();i++) v[i]=vv[i];
	delete [] vv;
#endif
}

void tpem_bcast(int *dS,TpemOptions const &tpemOptions)
{
#ifdef USE_MPI
	MPI_Bcast(dS,1,MPI_INT,0,tpemOptions.mpiCommTpem);       
	MPI_Barrier(tpemOptions.mpiCommTpem);
#endif
}

void tpem_bcast(double *dS,TpemOptions const &tpemOptions)
{
#ifdef USE_MPI
	MPI_Bcast(dS,1,MPI_DOUBLE,0,tpemOptions.mpiCommTpem);       
	MPI_Barrier(tpemOptions.mpiCommTpem);
#endif
}


double factor_alpha(int m)
{
	if (m==0) return 1;
	return 2;
}

double my_f (double x, void * p) {
   struct my_f_params * params 
     = (struct my_f_params *)p;

   /* return  params->funk(params->m,x);*/
   double tmp;
   double tmp2= (double)1.0/M_PI;
   int m=params->m;
   tmp = params->funk(x) *  tmp2 * factor_alpha(m) * tpem_chebyshev(m,x)/sqrt(1.0-x*x);
   return tmp;
}



void tpem_calculate_coeffs (unsigned int FKWCutOff,vector<double> &vobs, tpem_func_ptr funk)
{
	//int m;
	//struct my_f_params params;
	double *pts;
	//double result,result1,result2,abserr;
	//int npts,limit;
	//double epsabs,epsrel;
	
	
	
#ifndef NO_GSL	
	gsl_integration_workspace *workspace;
	gsl_function f;
	pts = new double[2];
	npts = 2;
	
	pts[0]= -1.0;
	pts[1] = 1.0;
	epsabs=1e-9;
	epsrel=1e-9;
	limit = 1e6;
	workspace =  gsl_integration_workspace_alloc (limit+2);

	result=result1=result2=abserr=0;

	params.funk = funk;
	f.function= &my_f;
	f.params = &params;
	for (m=0;m<FKWCutOff;m++) {
		params.m = m;
		gsl_integration_qagp (&f,pts,npts,epsabs,epsrel,limit,workspace,&result,&abserr);
		vobs[m] = result;
		
	}
	gsl_integration_workspace_free (workspace);
#else
	cerr<<"GSL is required to compute tpem_calculate_coeffs\n";
	cerr<<"Use tpem_calculate_coeffs_alt instead\n";
	cerr<<"At this point "<<__FILE__<<" "<<__LINE__<<endl;
	exit(1);
#endif
	delete [] pts;
}



void tpem_calculate_coeffs_alt (unsigned int FKWCutOff, vector<double> &vobs, tpem_func_ptr funk)
{
  unsigned int m,k;
  
  double *f;
  double x, sum;
  
  f = new double[FKWCutOff];
  
  for(m = 0; m<FKWCutOff; m++ ) 
    {
      x = cos(M_PI*(m+0.5)/FKWCutOff);
      f[m] = (*funk)(x);
    }
  
  for(m = 0; m<FKWCutOff; m++) 
    {
      sum = 0.0;
      for(k = 0; k<FKWCutOff; k++) 
	sum += f[k]*cos(M_PI * m*(k+0.5)/FKWCutOff);
      vobs[m] = 2.0 * sum/FKWCutOff;
    }
  vobs[0]/=2.0;
  delete [] f;
}

void tpem_calculate_coeffs (vector<double> &vobs, tpem_func_ptr funk,
	TpemOptions const &tpemOptions)
{
	if (tpemOptions.coeffs=="accurate" || tpemOptions.coeffs=="gsl") {
		tpem_calculate_coeffs(tpemOptions.cutoff,vobs,funk);
	} else {
		tpem_calculate_coeffs_alt(tpemOptions.cutoff,vobs,funk);
	}
}

double tpem_chebyshev(int m,double x) {
	double tmp;
	int p;
	if (m==0) return 1;
	if (m==1) return x;
	
	if ((m%2)==0) {
		p=m/2;
		tmp=tpem_chebyshev(p,x);
		return (2*tmp*tmp-1);
	}
	else {
		p=(m-1)/2;
		return (2*tpem_chebyshev(p,x)*tpem_chebyshev(p+1,x)-x);
	}
}	

void my_handler (const char * reason, const char * file, int line, int gsl_errno)
{
	cerr<<"GSL error handler called with reason="<<reason<<" file="<<file<<" line="<<line;
	cerr<<" gsl_errno="<<gsl_errno<<endl;
}


void tpem_off_diagonal_element_tpem (tpem_sparse *matrix, vector<vector<tpem_t> > &moment,
                unsigned int ket, TpemOptions const &tpemOptions)
{       /*
         *      on return:
         *      moment[bra][m] = <bra|Tm(X)|ket>
         *       
         */
	unsigned int n = tpemOptions.cutoff;
        tpem_subspace *info;
        vector<tpem_t>  tmp(matrix->rank), jm0(matrix->rank), jm1(matrix->rank);
        unsigned int  i, m, *p;
        tpem_t  keep;

        info = new_tpem_subspace (matrix->rank);
        
        for (m = 0; m < n; m++)
                for (i = 0; i < matrix->rank; i++)
                        moment[i][m] = 0.0;

        for (i = 0; i < matrix->rank; i++)
                tmp[i] = jm0[i] = jm1[i] = 0.0;
        jm0[ket] = 1;        /* set |j,0> */

        tpem_subspace_reset (info);
        tpem_subspace_push (info, ket);

        /* calculate |j,1> = X|j,0> */
	 tpem_sparse_product_tpem (matrix, info, jm1, jm0, tpemOptions,tpemOptions.epsProd);

        moment[ket][0] = jm0[ket];
        for (p = info->stack; p < info->top; p++)
                moment[*p][1] = jm1[*p];

        /* calculate |j,m> = 2X|j,m-1> - |j,m-2>
         *
         * begin (m=2) pass     jm0 = |j,0>     jm1 = |j,1>
         * end   (m=2) pass     jm0 = |j,1>     jm1 = |j,2>
         * ...
         * begin (m=k) pass     jm0 = |j,k-2>   jm1 = |j,k-1>
         * end   (m=k) pass     jm0 = |j,k-1>   jm1 = |j,k>
         */
        for (m = 2; m < n; m++) {
                /* calculate |tmp> = X|jm1> */
                tpem_sparse_product_tpem (matrix, info, tmp, jm1, tpemOptions,tpemOptions.epsProd);
                for (p = info->stack; p < info->top; p++) {
 			i = *p;
                        keep = tmp[i] + tmp[i] - jm0[i];
                        moment[i][m] = keep;
                        jm0[i] = jm1[i];
                        jm1[i] = keep;
                }
        }
	tpem_subspace_free(info);
}



void tpem_clean ()
{
	
}

void tpem_finalize()
{
	tpem_clean();
#ifdef USE_MPI
	int rank=0;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Finalize();
#endif
}


