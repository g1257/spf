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




#include "lanczosplus.h"

using namespace std;

void MatrixGet(MatType **matrix,int n,int pos,MatType &v)
{
    int x,y;
    y=int(n/pos);
    x=n % pos;
    v=matrix[x][y];
}

int SparseFromFile(SparseMatrix &mat, char *filename)
{
  int i, j, k;
  double eps = 1.e-15;
  MatType  Aij;
  int tmp;
  
  ifstream fin(filename);
  if (!fin || fin.bad()) {
    cerr<<"new_sparse: can't open file "<<filename<<endl;
    return -1;
  }

    fin>>tmp;
    mat.setSize(tmp);
    mat.initRow(mat.getSize()+1);

  k = 0;
  for (i = 0; i < mat.getSize(); i++) {
    mat.setRow(i,k);
    for (j = 0; j < mat.getSize(); j++) {
        fin>>Aij;
          if (abs (Aij) > eps)
        k++;
    }
  }
  mat.setRow(mat.getSize(),k); /* # of non-zero elements */

  mat.initCol(k);
  mat.initValues(k);
    
    fin.close();
    fin.open(filename);
    if (!fin || fin.bad()) {
    cerr<<"new_sparse: can't open file "<<filename<<endl;
    return -1;
  }

    fin>>tmp;
    mat.setSize(tmp);

  k = 0;
  for (i = 0; i < mat.getSize(); i++)
    for (j = 0; j < mat.getSize(); j++) {
      fin>>Aij;
      if (abs (Aij) > eps) {
    mat.setCol(k,j);
    mat.setValues(k,Aij);
    k++;
      }
    }

  fin.close();

  return mat.getSize();
}

int SparseFromMemory(SparseMatrix &mat, MatType **matrix,int n)
{
  int i, j, k,pos;
  double eps = 1.e-15;
  MatType  Aij;

    mat.setSize(n);
    mat.initRow(mat.getSize()+1);

  k = 0;
    pos=0;
  for (i = 0; i < mat.getSize(); i++) {
    mat.setRow(i,k);
    for (j = 0; j < mat.getSize(); j++) {
        MatrixGet(matrix,n,pos,Aij);
        pos++;
          if (abs (Aij) > eps)
        k++;
    }
  }
  mat.setRow(mat.getSize(),k); /* # of non-zero elements */

  mat.initCol(k);
  mat.initValues(k);
    
  pos=0;
  k = 0;
  for (i = 0; i < mat.getSize(); i++)
    for (j = 0; j < mat.getSize(); j++) {
      MatrixGet(matrix,n,pos,Aij);
      pos++;
      if (abs (Aij) > eps) {
    mat.setCol(k,j);
    mat.setValues(k,Aij);
    k++;
      }
    }

  return mat.getSize();
}





double ground (int flag, int n, vector<double> const &a, vector<double> const &b, 
vector<double> &gs)
{ /* adapted from tql1, handbook for automatic computation */
  // Set flag = 1 for the ground state
  // Set flag = 0 for the maximum eigenvalue

  int i, k, l, m;
  double c, dd, f, g, h, p, r, s, *d, *e, *v = 0, *vki;

  if (gs.size()>0) {
    v  = new double[n*n];
    for (k=0;k<n*n;k++) v[k]=0.0;
    for (k = 0, vki = v; k < n; k++, vki += (n + 1))
      (*vki) = 1.0;
  }

  d = new double[n]; 
  e = new double[n];

  for (i = 0; i < n; i++) {
    d[i] = a[i];
    e[i] = b[i];
  }

  for (l = 0; l < n; l++) {
    do {
      for (m = l; m < n - 1; m++) {
        dd = fabs (d[m]) + fabs (d[m + 1]);
        if ((fabs (e[m]) + dd) == dd) break;
      }
      if (m != l) {
        g = (d[l + 1] - d[l]) / (2.0 * e[l]);
        r = sqrt (g * g + 1.0);
        g = d[m] - d[l] + e[l] / (g + (g >= 0 ? fabs (r) : -fabs (r)));
        for (i = m - 1, s = c = 1.0, p = 0.0; i >= l; i--) {
          f = s * e[i];
          h = c * e[i];
          e[i + 1] = (r = sqrt (f * f + g * g));
          if (r == 0.0) {
            d[i + 1] -= p;
            e[m] = 0.0; 
            break;
          }
          s = f / r;
          c = g / r;
          g = d[i + 1] - p;
          r = (d[i] - g) * s + 2.0 * c * h;
          d[i + 1] = g + (p = s * r);
          g = c * r - h;
      if (gs.size()>0)
        for (k = 0, vki = v + i; k < n; k++, vki += n) {
          f = vki[1];
          vki[1] = s * vki[0] + c * f;
          vki[0] = c * vki[0] - s * f;
        }
        }
        if (r == 0.0 && i >= l) continue;
        d[l] -= p;
        e[l] = g;
        e[m] = 0.0;
      }
    } while (m != l);
  }

  if (flag)   // Search for minimal eigenvalue
    {
      for (i = 1, s = d[l = 0]; i < n; i++)
	if (d[i] < s) s = d[l = i];
    }
  else        // Search for maximal eigenvalue
    {
      for (i = 1, s = d[l = 0]; i < n; i++)
	if (d[i] > s) s = d[l = i];
    }

  if (gs.size()>0) {
    for (k = 0, vki = v + l; k < n; k++, vki += n)
      gs[k] = (*vki);
    free (v);
  }

    delete [] d;
    delete [] e;

  return s;
}

extern void vectorPrint(vector<double> const  &v,char *name,std::ostream &s);
extern void vectorPrint(vector<MatType> const &v,char *name,std::ostream &s);

int lanczos (int flag, SparseMatrix const &mat, int max_nstep, double &gsEnergy,
vector<MatType> &z)
{ /*
   *     In each step of the Lanczos algorithm the values of a[]
   *     and b[] are computed.
   *     then a tridiagonal matrix T[j] is formed from the matrix
   *     T[j-1] as
   *
   *            | a[0]  b[0]                            |
   *            | b[0]  a[1]  b[1]                      |
   *     T(j) = |       b[1]    .     .                 |
   *            |               .     .  a[j-2]  b[j-2] |
   *            |               .     .  b[j-2]  a[j-1] |
   */

  // Set flag = 1 for the ground state
  // Set flag = 0 for the maximum eigenvalue


  int i, j;
  long int seed;
  double s, enew, eold, atmp, btmp, ctmp, eps = 0.001;
  vector<double> a,b,c;
  vector<MatType> x,y;
  MatType tmp,tmp2;
    
	vector<double> nullVector;
	if (nullVector.size()>0) {
		std::cerr<<"nullVector has size greater than zero\n";
		exit(1);
	}
	
  seed= rand();    

  for (i = 0; i < max_nstep; i++) {
    a.push_back(0.0);
    b.push_back(0.0);
    c.push_back(0.0);
  }

  srand48 (seed);
  s = 0.0;
  z.clear();
  for (i = 0; i < mat.getSize(); i++) {
    tmp2 = MatType(0,0);
    z.push_back(tmp2);	 
    tmp2 = MatType(0,0);
    x.push_back(tmp2);
    tmp2 = MatType(drand48 () - 0.5, drand48 () - 0.5);
    y.push_back(tmp2);
    s += real (y[i] * conj (y[i]));
  }

  s = 1.0 / sqrt (s);
  for (i = 0; i < mat.getSize(); i++)
    y[i] *= s;

  if (max_nstep > mat.getSize())
    max_nstep = mat.getSize();

  eold = 100.;
   //cout<<"iter\t< 0 | H | 0 >\n";
  for (j = 0; j < max_nstep; j++) {
    mat.sparse_mult (x, y);

    atmp = 0.0;
    for (i = 0; i < mat.getSize(); i++)
      atmp += real (y[i] * conj (x[i]));

    btmp = 0.0;
    for (i = 0; i < mat.getSize(); i++) {
      x[i] -= atmp * y[i];
      btmp += real (x[i] * conj (x[i]));
    }

    a[j] = atmp;
    b[j] = (btmp = sqrt (btmp));
    enew = ground (flag,j+1, a, b,nullVector);
    //        cout<<j<<":\t"<<enew<<" "<<nullVector.size()<<endl;

    if (fabs ((enew - eold)/enew) < eps)
      break;

    eold = enew;
    for (i = 0; i < mat.getSize(); i++) {
      tmp = y[i];
      y[i] = x[i] / btmp;
      x[i] = -btmp * tmp;
    }
  }

  if (j < max_nstep)
    max_nstep = j + 1;

	
  gsEnergy = ground (flag,max_nstep, a, b, c);

  srand48 (seed);
  for (i = 0; i < mat.getSize(); i++) {
    x[i] = MatType(0,0);    
    y[i] = s * MatType((drand48 () - 0.5), (drand48 () - 0.5)) ;
   
    z[i] = MatType(0,0);
     
  }

  for (j = 0; j < max_nstep; j++) {
 
    mat.sparse_mult (x, y);
 
    
    atmp = a[j];
    btmp = b[j];
    ctmp = c[j];
 
    for (i = 0; i < mat.getSize(); i++) {
      z[i] += ctmp * y[i];
      x[i] -= atmp * y[i];
      tmp = y[i];
      y[i] = x[i] / btmp;
      x[i] = -btmp * tmp;
     
    }
  }

  return 0;
}

void initMatrix(MatType **matrix,int n)
{
    double temp,max_value=10;
    int i,j;
    
    srand(13238);
    for (i=0;i<n;i++) {
        for (j=0;j<i;j++) {
            matrix[i][j]=matrix[j][i]=0;
        }
    }
    for (i=0;i<n;i++) {
        for (j=0;j<i;j++) {
            if ((double) rand()/RAND_MAX < 0.05 ) {
                temp=max_value*(double)rand()/RAND_MAX;
                matrix[i][j]=matrix[j][i]=temp;
            }
        }
    }
}



void check1(SparseMatrix const &matrix,int n,vector<MatType> const &gsVector)
{
    int i;
    MatType s=0;
    vector<MatType> x;
    
    
    for (i=0;i<n;i++) x.push_back(0.0);
    matrix.sparse_mult(x,gsVector);
    for (i=0;i<n;i++) {
        s+=x[i]*conj(gsVector[i]);
    }
    cout<<"Check1: "<<s<<endl;
}
    
