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




#include "common.h"
#include "io.h"
#include "Matrix.h"

#ifdef USE_MPI
extern
void MPI_VecReduce(vector<tpem_t> const &v1,vector<tpem_t> &v2,int n,MPI_Datatype type1,MPI_Op op1,int rank,MPI_Comm comm);
#endif


extern double Density(double x);

////////%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


#ifdef MODEL_KONDO_DMS_MANYBANDS
double calcNumber(DynVars const &dynVars,Geometry const &geometry,Parameters const &ether,Aux &aux,int band)
{
	int lambda,i,N=geometry.volume();
	double s=0;
	MatType temp;
	
	diag(aux.eigOneBand,geometry,dynVars,ether,aux,'V');
	for (lambda=0;lambda<ether.hilbertSize;lambda++) {
			temp=0;
			for (i=0;i<N;i++) {
				temp+=conj(aux.matrix(i+N*band,lambda))*aux.matrix(i+N*band,lambda);
			}
			s+=real(temp)*fermi(ether.beta*(aux.eigOneBand[lambda]-aux.varMu)); 
	}
	return s;

}
#endif



#ifdef  MODEL_KONDO_DMS_CUBIC
void accNOfOmega(Geometry const &geometry,DynVars const &dynVars, Parameters const &ether,Aux &aux)
{
	int lambda,i,N=geometry.volume(),spin;
	double weight;
	MatType temp;
	double shiftmu=0.0;
	
	if (ether.isSet("shiftmu")) {
		shiftmu = aux.varMu;
	}
	diag(aux.eigOneBand,geometry,dynVars,ether,aux,'V');
	for (lambda=0;lambda<ether.hilbertSize;lambda++) {
		for (spin=0;spin<2;spin++) {
			temp=0;
			for (i=0;i<N;i++) {
				temp+=conj(aux.matrix(i+spin*N,lambda))*aux.matrix(i+spin*N,lambda);
			}
			weight=real(temp);
			aux.Nw[0].add(aux.eigOneBand[lambda]-shiftmu,weight); // dos a
		}
	}
	for (lambda=0;lambda<ether.hilbertSize;lambda++) {
		for (spin=0;spin<2;spin++) {
			temp=0;
			for (i=0;i<N;i++) {
				temp+=conj(aux.matrix(i+spin*N+2*N,lambda))*aux.matrix(i+spin*N+2*N,lambda);
			}
			weight=real(temp);
			aux.Nw[1].add(aux.eigOneBand[lambda]-shiftmu,weight); //dos b
		}
	}
	
}
#endif

#ifdef  MODEL_KONDO_DMS_ZINCBLENDE
void accNOfOmega(Geometry const &geometry,DynVars const &dynVars, Parameters const &ether,Aux &aux)
{
	int lambda,i,N=geometry.volume(),spin;
	double weight;
	MatType temp;
	double shiftmu=0.0;
	
	if (ether.isSet("shiftmu")) {
		shiftmu = aux.varMu;
	}
	diag(aux.eigOneBand,geometry,dynVars,ether,aux,'V');
	for (orbital=0;orbital<4;orbital++) {
		for (lambda=0;lambda<ether.hilbertSize;lambda++) {	
			temp=0;
			for (i=0;i<N;i++) {
				temp+=conj(aux.matrix(i+orbital*N,lambda))*aux.matrix(i+orbital*N,lambda);
			}
			weight=real(temp);
		}
		aux.Nw[orbital].add(aux.eigOneBand[lambda]-shiftmu,weight); // dos for orbital
	}
	
}
#endif

#ifdef  MODEL_KONDO_DMS_FCC


void accNOfOmega(Geometry const &geometry,DynVars const &dynVars, Parameters const &ether,Aux &aux)
{
	int lambda,i,N=geometry.volume(),spin;
	double weight;
	MatType temp;
	double shiftmu=0.0;
	int orbital;
	
	if (ether.isSet("shiftmu")) {
		shiftmu = aux.varMu;
	}
	diag(aux.eigOneBand,geometry,dynVars,ether,aux,'V');
	for (orbital=0;orbital<6;orbital++) {
		for (lambda=0;lambda<ether.hilbertSize;lambda++) {
			temp=0;
			for (i=0;i<N;i++) {
				temp+=conj(aux.matrix(i+orbital*N,lambda))*aux.matrix(i+orbital*N,lambda);
			}
			weight =real(temp);
			aux.Nw[orbital].add(aux.eigOneBand[lambda]-shiftmu,weight); // dos for orbital
		}
		
	}
	
}
#endif


#ifndef MODEL_KONDO_DMS_MANYBANDS
void accNOfOmega(Geometry const &geometry,DynVars const &dynVars, Parameters const &ether,Aux &aux)
{
	int lambda,i,N=geometry.volume();
	double weight;
	MatType temp;
	double shiftmu=0.0;
	
	if (ether.isSet("shiftmu")) {
		shiftmu = aux.varMu;
	}
	//diag(aux.eigOneBand,geometry,dynVars,ether,aux,'V');
	for (lambda=0;lambda<ether.hilbertSize;lambda++) {
		temp=0;
		for (i=0;i<ether.hilbertSize;i++) {
			temp+=conj(aux.matrix(i,lambda))*aux.matrix(i,lambda);
		}
		weight=real(temp);
		aux.Nw[0].add(aux.eigOneBand[lambda]-shiftmu,weight);
	}
}

#endif
//////////%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

void accAkw(Geometry const &geometry,DynVars const &dynVars, Parameters const &ether,Aux &aux)
{
	int r,lambda,i,j;
	int n=ether.linSize,int_dof=ether.hilbertSize/ether.linSize;
	MatType temp;
	
	
	//diag(aux.eigOneBand,geometry,dynVars,ether,aux,'V');
#ifdef MODEL_KONDO_FINITE
	int spin;
	for (r=0;r<n;r++) {
		for (lambda=0;lambda<ether.hilbertSize;lambda++) {
			temp=0.0;
			for (i=0;i<n;i++) {
				j=geometry.add(i,r);
				for (spin=0;spin<int_dof;spin++) 
					temp += conj(aux.matrix(i+spin*n,lambda))*aux.matrix(j+spin*n,lambda);
			}
			
			aux.Arw[r].add(aux.eigOneBand[lambda],real(temp));
			aux.ArwC[r].add(aux.eigOneBand[lambda],imag(temp));
		}
	}

#else
#ifdef MODEL_KONDO_INF_ONEBAND	
	for (r=0;r<n;r++) {
		for (lambda=0;lambda<ether.hilbertSize;lambda++) {
			temp=0.0;
			for (i=0;i<n;i++) {
				j=geometry.add(i,r);
				temp += conj(aux.matrix(i,lambda))*aux.matrix(j,lambda);
			}
			
			aux.Arw[r].add(aux.eigOneBand[lambda],real(temp));
			aux.ArwC[r].add(aux.eigOneBand[lambda],imag(temp)); // -aux.varMu optional
		}
	}
#else
#ifdef MODEL_KONDO_DMS_FCC
	// here A(r+gamma*N,omega) contains A(r,omega)_\gamma
	int gamma;
	for (r=0;r<n;r++) {
		for (lambda=0;lambda<ether.hilbertSize;lambda++) {
			for (gamma=0;gamma<6;gamma++) {
				temp=0.0;
				for (i=0;i<n;i++) {
					j=geometry.add(i,r);
					temp += conj(aux.matrix(i+gamma*n,lambda))*aux.matrix(j+gamma*n,lambda);
					aux.Arw[r+gamma*n].add(aux.eigOneBand[lambda],real(temp));
					aux.ArwC[r+gamma*n].add(aux.eigOneBand[lambda],imag(temp));
				}
			}
		}
	}
#else

	cerr<<"Calculation of A(k,omega) is not implemented for this model (sorry) "<<__LINE__<<endl;
#endif
#endif
#endif
}


		
void accLcd(Geometry const &geometry,DynVars const &dynVars,Parameters const &ether,Aux &aux)
{
	int i,lambda,alpha,int_dof=ether.hilbertSize/ether.linSize;
	MatType tmp;
	double s;
	

	for (i=0;i<ether.linSize;i++) {
		s=0.0;
		for (lambda=0;lambda<ether.hilbertSize;lambda++) {
			for (alpha=0;alpha<int_dof;alpha++) {
				tmp = conj(aux.matrix(i+alpha*ether.linSize,lambda))*aux.matrix(i+alpha*ether.linSize,lambda);
				s+= real(tmp)*fermi(ether.beta*(aux.eigOneBand[lambda]-aux.varMu));
			}
		}
		if (ether.isSet("savelcd")) {
			aux.lcd[i]=s;
		} else {
			//aux.lcd[i] += s;
		}
	}	
}

double calcKinetic(DynVars const &dynVars,Geometry const &geometry,Parameters const &ether,Aux &aux)
{
	int i,j,k,lambda,border;
	double ret=0.0,tmp2;
	MatType tmp,hopping;
		
	//diag(aux.eigOneBand,geometry,dynVars,ether,aux,'V');
	
#ifdef MODEL_KONDO_FINITE
	int spin;
	for (lambda=0;lambda<ether.hilbertSize;lambda++) {
		tmp2=0.0;
	for (i=0;i<ether.linSize;i++) {
		for (k=0;k<geometry.z(i);k++) {
			j=geometry.neighbor(i,k);
			border=geometry.borderId(i,k);
			if (border>=0) hopping = ether.hoppings[border];
			else hopping = -1.0;
			for (spin=0;spin<2;spin++) {
				tmp = conj(aux.matrix(i+spin*ether.linSize,lambda))*aux.matrix(j+spin*ether.linSize,lambda);
				tmp2 += real(hopping)*real(tmp);
			}
		}
	}
		ret += tmp2 * fermi(ether.beta*(aux.eigOneBand[lambda]-aux.varMu));
	}
	
#else
#ifdef MODEL_KONDO_INF_ONEBAND
	for (lambda=0;lambda<ether.hilbertSize;lambda++) {
		tmp2=0.0;
		for (i=0;i<ether.linSize;i++) {
			for (k=0;k<geometry.z(i);k++) {
				j=geometry.neighbor(i,k);
				border=geometry.borderId(i,k);
				if (border>=0) hopping = ether.hoppings[border];
				else hopping = -1.0;
				tmp = conj(aux.matrix(i,lambda))*aux.matrix(j,lambda);
				tmp2 += real(hopping)*real(tmp);
			}
		}
		ret += tmp2 * fermi(ether.beta*(aux.eigOneBand[lambda]-aux.varMu));
	}
#else
	ret= -1; // not implemented
#endif	
#endif
	return ret;
}

void vectorAcc(vector<double> &sm,vector<tpem_t> const &moment,int size)
{
	int i;
	for (i=0;i<size;i++) {
		sm[i] += real(moment[i]);
	}
}

double measure_kinetic(Geometry const &geometry,Parameters const &ether,Aux &aux,TpemOptions const &tpemOptions)
{
	int cutoff = ether.tpem_cutoff;
	static int firstcall=1;
	static vector<double>	coeffs(cutoff);
	static vector<vector<tpem_t> > moment;
	static vector<tpem_t> tmpVector(cutoff);
	int i,j,k,border;
	MatType hopping;
	double beta=ether.beta;
	vector<double> save_moment(cutoff);
	double ret=0.0;
	
	for (i=0;i<cutoff;i++) save_moment[i]=0.0;
	
	if (firstcall || ether.carriers>0) {
		firstcall=0;
		tmpValues(aux.varTpem_a,aux.varTpem_b,aux.varMu,beta,0);
		tpem_calculate_coeffs (coeffs, Density,tpemOptions);
		for (i=0;i<ether.hilbertSize;i++) moment.push_back(tmpVector);
	}
	tmpValues(aux.varTpem_a,aux.varTpem_b,aux.varMu,beta,0);
	if (ether.isSet("adjusttpembounds")) tpem_calculate_coeffs (coeffs, Density,tpemOptions);
	
	
#ifdef MODEL_KONDO_FINITE
	int spin;
	for (i=0;i<ether.linSize;i++) {
		for (k=0;k<geometry.z(i);k++) {
			j = geometry.neighbor(i,k);
			border=geometry.borderId(i,k);
			if (border>=0) hopping = ether.hoppings[border];
			else hopping = -1.0;
			for (spin=0;spin<2;spin++) {
		tpem_off_diagonal_element_tpem (aux.sparseMatrix[0], moment,i+spin*ether.linSize,tpemOptions);
		ret += real(hopping)*tpem_expansion_complex (moment[j+spin*ether.linSize], coeffs);
			if (ether.isSet("dumpkinetic")) vectorAcc(save_moment,moment[j+spin*ether.linSize],cutoff);
			}
		}
	}
	if (ether.isSet("dumpkinetic")) {
		string s = ether.rootname;
		s =s + ".kin"; 
		std::ofstream fout(s.c_str());
		for (i=0;i<cutoff;i++) {
			fout<<i<<" "<<coeffs[i]<<" "<<save_moment[i]<<endl;
		}
		fout.close();
	}
	
#else
#ifdef MODEL_KONDO_INF_ONEBAND
	for (i=0;i<ether.linSize;i++) {
		for (k=0;k<geometry.z(i);k++) {
			j = geometry.neighbor(i,k);
			border=geometry.borderId(i,k);
			if (border>=0) hopping = ether.hoppings[border];
			else hopping = -1.0;
			tpem_off_diagonal_element_tpem (aux.sparseMatrix[0],moment,i,tpemOptions);
			ret += real(hopping)*tpem_expansion_complex (moment[j], coeffs);
			
		}
	}
#else	
	ret= -1; // not implemented
#endif
#endif
	return ret;
}	


void kTpemMomentsOptical(Geometry const &geometry, Parameters const &ether,Aux &aux,TpemOptions const &tpemOptions)
{
	int cutoff=ether.tpem_cutoff;
	static int firstcall=1;
	static vector<vector<tpem_t> > moment,moment2;
	static vector<tpem_t> tmpVector(cutoff+5);
	int m,l,i,j,i2,j2,sign;
	int n=ether.linSize;
	vector<int> nx(2*n);
	int tmp;
	tpem_t val;
	unsigned int k,activeProc,part;
	double buf0,buf1;
	
	
	if (firstcall) {
		firstcall=0;
		for (i=0;i<ether.hilbertSize;i++) moment.push_back(tmpVector);
		for (i=0;i<ether.hilbertSize;i++) moment2.push_back(tmpVector);
		
		
		for (i=0;i<n;i++) {
			nx[i] = geometry.neighbor(i,0); // i+x
			nx[i+n] = geometry.neighbor(i,1); // i-x
		}
		
		for (m=0;m<cutoff;m++)  for (l=0;l<cutoff;l++) aux.opticalMoments.push_back(0.0);
		
	}
	
	
#ifdef MODEL_KONDO_FINITE
	part=ether.hilbertSize/ether.mpiNop1;
        if (ether.hilbertSize % ether.mpiNop1>0) part++;
        activeProc=ether.hilbertSize/part;
        if (ether.hilbertSize % part>0) activeProc++;
	
	for (k=0;k<part;k++) {
		i=k+ether.mpiRank * part;
		tpem_off_diagonal_element_tpem (aux.sparseMatrix[0], moment,i,tpemOptions);
		for (i2=0;i2<2;i2++) {
			if (i<n) tmp=nx[i+i2*n];
			else tmp=nx[i-n+i2*n]+n;
			tpem_off_diagonal_element_tpem (aux.sparseMatrix[0], moment2,tmp,tpemOptions);
			for (l=0;l<cutoff;l++) {
				for (m=0;m<cutoff;m++) {
					for (j=0;j<ether.hilbertSize;j++) {
						for (j2=0;j2<2;j2++) {
							if (j<n) tmp=nx[j+j2*n];
							else tmp=nx[j-n+j2*n]+n;
							sign = (2*i2 -1)*(2*j2-1);
							val = moment[j][l] * moment2[tmp][m];
							aux.opticalMoments[l+cutoff*m] += double(sign)*val;
						}
					}
				}
			}
		}
	}

#ifdef USE_MPI	
        MPI_Status status;
        MPI_VecReduce(aux.opticalMoments,aux.opticalMoments,cutoff*cutoff,MPI_DOUBLE,MPI_SUM,0,tpemOptions.mpiCommTpem);

#else
	if (ether.mpiRank!=0 || ether.mpiNop1!=1) {
		cerr<<"rank should be zero and nop should be 1 for serial version but are ";
		cerr<<ether.mpiRank<<" and "<<ether.mpiNop1<<endl;
		exit(1);
	}																			
#endif //use_mpi
#else
#ifdef MODEL_KONDO_INF_ONEBAND
	// parallelize this code (FIXME)
	for (i=0;i<ether.hilbertSize;i++) {
		tpem_off_diagonal_element_tpem (aux.sparseMatrix[0], moment,i,tpemOptions);
		for (i2=0;i2<2;i2++) {
			tpem_off_diagonal_element_tpem (aux.sparseMatrix[0],moment2,nx[i+n*i2],tpemOptions);
			for (l=0;l<cutoff;l++) {
				for (m=0;m<cutoff;m++) {
					for (j=0;j<ether.hilbertSize;j++) {
						for (j2=0;j2<2;j2++) {
							sign = (2*i2 -1)*(2*j2-1);
							val = moment[j][l] * moment2[nx[j+n*j2]][m];
							aux.opticalMoments[l+cutoff*m] =
aux.opticalMoments[l+cutoff*m]+double(sign)*val;
						}
					}
				}
			}
		}
	}
	
#else
	cerr<<"Calculation of sigma(k,omega) is not implemented for this model (sorry) "<<__LINE__<<endl;
#endif
#endif
}

#ifdef MODEL_KONDO_DMS_FCC
void accOptical(Geometry const &geometry,DynVars const &dynVars,Parameters const &ether,Aux &aux)
{
	int lambda2,lambda,i,j,spin,spin2;
	int n=ether.linSize,int_dof=ether.hilbertSize/ether.linSize;
	MatType temp,thop;
	double e1,e2,beta=ether.beta;
	//int aindex = 0*(ether.linSize/4); //(0,0,0)
	//int aindex = 1*(ether.linSize/4); //(0,a/2,a/2)
	//int aindex = 2*(ether.linSize/4); //(a/2,0,a/2)
	int aindex = 3*(ether.linSize/4); //(a/2,a/2,0)
	
	int dir;
	
	//diag(aux.eigOneBand,geometry,dynVars,ether,aux,'V');
	// some checking
	if (geometry.latticeName()!="fcc" || ether.linSize % 4!=0) {
		cerr<<"INTERNAL AT "<<__FILE__<<" "<<__LINE__<<endl;
		exit(1);
	}
	
	for (lambda=0;lambda<ether.hilbertSize;lambda++) {
		e1 = aux.eigOneBand[lambda]-aux.varMu;
		for (lambda2=0;lambda2<ether.hilbertSize;lambda2++) {			
			e2 = aux.eigOneBand[lambda2]-aux.varMu;
			if (e1-e2<1e-3) continue; 
			temp = MatType(0,0);
			for (i=0;i<n;i++) {
				j = geometry.add(i,aindex); // j = i+avector
				dir = geometry.scalarDirection(i,j,0); // third entry does not matter here
				for (spin=0;spin<int_dof;spin++) {
					for (spin2=0;spin2<int_dof;spin2++) {
						thop = ether.bandHoppings[spin+spin2*4+16*dir];
						temp += thop * conj(aux.matrix(i+spin*n,lambda))*aux.matrix(j+spin2*n,lambda2);
						temp -= thop * conj(aux.matrix(j+spin*n,lambda)) * aux.matrix(i+spin2*n,lambda2);
					}
				}
			}
			temp = real(temp)*real(temp)+imag(temp)*imag(temp);
			
			temp = temp *(fermi(beta*e2)-fermi(beta*e1))/(e1-e2);
			//temp = temp * fermi(beta*e2) * fermi(-beta*e1);
			//temp = temp * (1.0 - exp(-beta*(e1-e2)))/(e1-e2);
			aux.Sigma.add(e1-e2,real(temp));
		}
	}
}	
#else
void accOptical(Geometry const &geometry,DynVars const &dynVars,Parameters const &ether,Aux &aux)
{
	
	int lambda2,lambda,i,j,spin;
	int n=ether.linSize,int_dof=ether.hilbertSize/ether.linSize;
	MatType temp;
	double e1,e2,beta=ether.beta;
	
	//diag(aux.eigOneBand,geometry,dynVars,ether,aux,'V');
	
	
	for (lambda=0;lambda<ether.hilbertSize;lambda++) {
		e1 = aux.eigOneBand[lambda]-aux.varMu;
		for (lambda2=0;lambda2<ether.hilbertSize;lambda2++) {			
			e2 = aux.eigOneBand[lambda2]-aux.varMu;
			if (e1-e2<1e-3) continue; 
			temp = MatType(0,0);
			for (i=0;i<n;i++) {
				j = geometry.neighbor(i,0); // j = i+x
				for (spin=0;spin<int_dof;spin++) {
					temp += conj(aux.matrix(i+spin*n,lambda))*aux.matrix(j+spin*n,lambda2);
					temp -= conj(aux.matrix(j+spin*n,lambda)) * aux.matrix(i+spin*n,lambda2);
				}
			}
			temp = real(temp)*real(temp)+imag(temp)*imag(temp);
			
			temp = temp *(fermi(beta*e2)-fermi(beta*e1))/(e1-e2);
			//temp = temp * fermi(beta*e2) * fermi(-beta*e1);
			//temp = temp * (1.0 - exp(-beta*(e1-e2)))/(e1-e2);
			aux.Sigma.add(e1-e2,real(temp));
		}
	}
}
#endif

void kTpemMomentsCl(Geometry const &geometry, Parameters const &ether,Aux &aux,TpemOptions const &tpemOptions)
{
	
	int cutoff=ether.tpem_cutoff;
	static int firstcall=1;
	static vector<vector<tpem_t>  > moment;
	vector<tpem_t> tmpVector(tpemOptions.cutoff);
	int r,lambda,i,j,spin,m;
	int n=ether.linSize,int_dof=ether.hilbertSize/ether.linSize;
	tpem_t tmp;
	unsigned int k,activeProc,part;
	double buf0,buf1;
	
	if (firstcall) {
		firstcall=0;
		for (i=0;i<ether.hilbertSize;i++) moment.push_back(tmpVector);
	}
	

#ifdef MODEL_KONDO_FINITE
	part=n/ether.mpiNop1;
        if (n % ether.mpiNop1>0) part++;
        activeProc=n/part;
        if (n % part>0) activeProc++;
	
	for (k=0;k<part;k++) {
		i=k + ether.mpiRank * part;
		for (spin=0;spin<int_dof;spin++) {					
			tpem_off_diagonal_element_tpem (aux.sparseMatrix[0], moment,i+spin*ether.linSize,
			tpemOptions);
			for (lambda=0;lambda<ether.hilbertSize;lambda++) {
				for (r=0;r<n;r++) {
					j=geometry.add(i,r);
					for (m=0;m<cutoff;m++) {
						//cerr<<"r="<<r<<" lambda="<<lambda<<" i="<<i<<" j="<<j<<" spin="<<spin<<" m="<<m<<endl;
						tmp = moment[j+spin*ether.linSize][m];
						aux.offdMoments[r+m*n] =aux.offdMoments[r+m*n] + tmp;
					}
				}
			}
		}
	}

#ifdef USE_MPI
	MPI_Status status;
        MPI_VecReduce(aux.offdMoments,aux.offdMoments,n*cutoff,MPI_DOUBLE,MPI_SUM,0,tpemOptions.mpiCommTpem);	
#else
	if (ether.mpiRank!=0 || ether.mpiNop1!=1) {
		cerr<<"rank should be zero and nop should be 1 for serial version but are "<<ether.mpiRank;
		cerr<<" and "<<ether.mpiNop1<<endl;
		exit(1);
	}																			
#endif //use_mpi			
#else
#ifdef MODEL_KONDO_INF_ONEBAND
	// fixme: parallelize this code (FIXME)
	for (i=0;i<n;i++) {				
		tpem_off_diagonal_element_tpem (aux.sparseMatrix[0],moment,i,tpemOptions);
		for (lambda=0;lambda<ether.hilbertSize;lambda++) {
			for (r=0;r<n;r++) {
				j=geometry.add(i,r);
				for (m=0;m<cutoff;m++) {
					//cerr<<"r="<<r<<" lambda="<<lambda<<" i="<<i<<" j="<<j<<" spin="<<spin<<" m="<<m<<endl;
					tmp = moment[j][m];
					aux.offdMoments[r+m*n] =aux.offdMoments[r+m*n] + tmp;
				}
			}
		}
	}
#else
#ifdef MODEL_KONDO_DMS_FCC
	part=n/ether.mpiNop1;
        if (n % ether.mpiNop1>0) part++;
        activeProc=n/part;
        if (n % part>0) activeProc++;
	
	for (k=0;k<part;k++) {
		i=k + ether.mpiRank * part;
		for (spin=0;spin<6;spin++) {	// "spin" is in reality internal degree of freedom				
			tpem_off_diagonal_element_tpem (aux.sparseMatrix[0], moment,i+spin*ether.linSize,
			tpemOptions);
			for (lambda=0;lambda<ether.hilbertSize;lambda++) {
				for (r=0;r<n;r++) {
					j=geometry.add(i,r);
					for (m=0;m<cutoff;m++) {
						aux.offdMoments[r+m*n] =aux.offdMoments[r+m*n]+moment[j+spin*ether.linSize][m];
					}
				}
			}
		}
	}

#ifdef USE_MPI
	MPI_Status status;
        MPI_VecReduce(aux.offdMoments,aux.offdMoments,n*cutoff,MPI_DOUBLE,MPI_SUM,0,tpemOptions.mpiCommTpem);	
#else
	if (ether.mpiRank!=0 || ether.mpiNop1!=1) {
		cerr<<"rank should be zero and nop should be 1 for serial version but are "<<ether.mpiRank;
		cerr<<" and "<<ether.mpiNop1<<endl;
		exit(1);
	}																			
#endif //use_mpi
#else
	cerr<<"Calculation of A(k,omega) is not implemented for this model (sorry) "<<__LINE__<<endl;
#endif
#endif
#endif
}

void accLdos(Geometry const &geometry,DynVars const &dynVars, Parameters const &ether,Aux &aux)
{
	int i,lambda,spin,n=geometry.volume();
	MatType temp;
	double shiftmu=0;
	
	
	if (ether.isSet("shiftmu")) {
		shiftmu = aux.varMu;
	}
	for (i=0;i<n;i++) {
		temp=0;
		for (lambda=0;lambda<ether.hilbertSize;lambda++) {
			for (spin=0;spin<2;spin++) {
				temp+=conj(aux.matrix(i+spin*n,lambda))*aux.matrix(i+spin*n,lambda);
			}
		}
		aux.Ldos[i].add(aux.eigOneBand[lambda]-shiftmu,real(temp));
	}
}

// verif. May 2006
MatType greenFunction(int index1,int index2, Parameters const &ether,Aux &aux)
{
	MatType tmp;
	int lambda;
	
	tmp =0;
	for (lambda=0;lambda<ether.hilbertSize;lambda++) {
		tmp += conj(aux.matrix(index1,lambda))*aux.matrix(index2,lambda)*fermi(-ether.beta*(aux.eigOneBand[lambda]-aux.varMu));
	}
	return tmp;
}


// New charge correlation function. Nov. 2005, verif. May 2006
void accChargeCorrelation(size_t gamma, size_t gamma2, 
                                Geometry const &geometry,DynVars const &dynVars, Parameters const &ether,Aux &aux)
{
        size_t x,w,v; //gamma is the orbital index
        int n = geometry.volume();
        MatType tmp;
        int int_dof =  int(ether.hilbertSize/n);
	double density = 0.5*0.75; //density per orbital per site, hence the 0.5 factor
		
        // Do for all x's
        for (x=0;x<n;x++) {     // Sum over all i's (or w's)
                tmp=0;
                for (w=0;w<n;w++) {
                        //for (gamma=0;gamma<int_dof;gamma++) {
                             //   for (gamma2=0;gamma2<int_dof;gamma2++) {
                                        v=geometry.add(x,w);
                                        if (x==0 && gamma==gamma2) {
                                                tmp += (1.0-greenFunction(w+gamma*n,w+gamma*n,ether,aux));
                                        } else {
						tmp +=(1.0-greenFunction(w+gamma*n,w+gamma*n,ether,aux))*
						(1.0-greenFunction(v+gamma2*n,v+gamma2*n,ether,aux));
                                                tmp -= greenFunction(v+gamma2*n,w+gamma*n,ether,aux)*
						greenFunction(w+gamma*n,v+gamma2*n,ether,aux);
                                      //  }
                               // }
                        }
                }
                aux.cco[x+gamma*n+gamma2*n*int_dof] += real(tmp)/n;
		aux.cco[x+gamma*n+gamma2*n*int_dof] -= density*density;
        }
}


// For MODEL_KONDO_INF_TWOBANDS
// T_{g1,g2} where g1 and g2 are the orbitals
MatType calcOrbitalT(int i,int g1, int g2,Parameters const &ether,Aux &aux)
{
	MatType scomplex=0.0;
	int lambda,volume = ether.linSize;
	
	for (lambda=0;lambda<2;lambda++) {
		scomplex+=conj(aux.matrix(i+g1*volume,lambda)) * aux.matrix(i+g2*volume,lambda)*
			fermi(ether.beta*(aux.eigOneBand[lambda]-aux.varMu));
	}
	return scomplex;
}
		
// For MODEL_KONDO_INF_TWOBANDS	
// Phys. Reports 344,1 Eq. (10) but there's no spin here
// dir=0,1,2 corresponding to x,y,z		
MatType calcOrbitalT(int i,int dir,Parameters const &ether,Aux &aux)
{
	switch (dir) {
		case 0: // T_i^x
			return 0.5*(calcOrbitalT(i,0,1,ether,aux)+calcOrbitalT(i,1,0,ether,aux));
			break;
		case 1: // T_i^y
			return 0.5*(calcOrbitalT(i,0,1,ether,aux)-calcOrbitalT(i,1,0,ether,aux));
			break;
		case 2: // T_i^z
			return 0.5*(calcOrbitalT(i,0,0,ether,aux)-calcOrbitalT(i,1,1,ether,aux));
			break;
	}
	if (ether.mpiRank==0) { // only master prints error message
		cerr<<"calcOrbitalT: unknown direction dir="<<dir<<endl;
		cerr<<"AT THIS POINT: "<<__FILE__<<" "<<__LINE__<<endl;
	}
	exit(1);
	return 0.0; // to avoid compiler warnings
}
			
// For MODEL_KONDO_INF_TWOBANDS	
// orbital angles as defined by H. Aliaga, Feb, 2006	
void accOrbitalAngles(Geometry const &geometry,Parameters const &ether,Aux &aux)
{
	int i, volume = ether.linSize;
	double tx,tz,psi;
	
	cerr<<"Accing orbital angles...\n";
	for (i=0;i<volume;i++) {
		tz=real(calcOrbitalT(i,0,ether,aux)); // "2" instead of "0"?
		tx=real(calcOrbitalT(i,2,ether,aux)); // "0" instead of "2"?
		cerr<<"CHECK "<<tx<<" "<<tz<<endl;
		if (fabs(tz)>1e-10) {
			psi = atan(tx/tz);
			psi += (1.0-mySign(tz))*M_PI*0.5;
		} else {
			psi = mySign(tx)*M_PI*0.5;
		}
		aux.orbitalAngles[i] += psi;
	}
}

extern double calcPhonon(int ind,DynVars const &dynVars,Geometry const &g,Parameters const &ether,int what);

MatType calcF(int i, int g1, int g2, DynVars const &dynVars, Geometry const &g,Parameters const &ether, Aux &aux)
{
	//tz = real(calcOrbitalT(i,0,ether,aux));
	//tx = real(calcOrbitalT(i,2,ether,aux));
	double q2 = calcPhonon(i,dynVars,g,ether,1);
	double q3 = calcPhonon(i,dynVars,g,ether,2);	
	double ksi = atan(q2/q3);
	
	if (g1 == 0  && g2 == 0) {
		return cos(ksi*0.5);
		//return cos(0.5*ksi)*cos(0.5*ksi);
	} else if (g1 == 0 && g2 == 1) {
		return sin(ksi*0.5);
		//return cos(0.5*ksi)*sin(0.5*ksi);
	} else  if (g1 == 1 && g2 == 0) {
		return -sin(0.5*ksi);
		//return -cos(0.5*ksi)*sin(0.5*ksi);
	} else if (g1 == 1 && g2 == 1) {
		return cos(0.5*ksi);
	} else {
		cerr<<"3rd orbital non-existent"<<endl;
		exit(1);
		// throw std::runtime_error("string"); // #include <stdexcept>
	}
}

void setupGreenFunction(psimag::Matrix<MatType> &G,Parameters const &ether,Aux &aux)
{
	for (size_t i=0;i<G.n_row();i++) for (size_t j=0;j<G.n_col();j++) G(i,j) = greenFunction(i,j,ether,aux);
}

MatType setupGtildeAux(size_t i,size_t gamma,size_t j,size_t gamma2,const psimag::Matrix<MatType>& G,
		DynVars const &dynVars,Geometry const &g,Parameters const &ether,Aux &aux)
{
	MatType sum = 0;
	size_t n = ether.linSize;
	size_t int_dof =  size_t(ether.hilbertSize/n);
	for (size_t q=0;q<int_dof;q++) {
		for (size_t q2=0;q2<int_dof;q2++) {
			sum += calcF(i,gamma,q,dynVars,g,ether,aux) * G(i+q*n,j+q2*n) * conj(calcF(j,gamma2,q2,dynVars,g,ether,aux));
			//sum += conj(calcF(i,q,gamma,ether,aux)) * G(i+q*n,j+q2*n) * calcF(j,q2,gamma2);
		}
	}
	return sum;
}

void setupGtilde(psimag::Matrix<MatType> &Gtilde,const psimag::Matrix<MatType>& G,DynVars const &dynVars,Geometry const &g,Parameters const &ether,Aux &aux)
{
	size_t n = ether.linSize;
	size_t int_dof =  size_t(ether.hilbertSize/n);
	
	for (size_t i=0;i<n;i++) {
		for (size_t j=0;j<n;j++) {
			 for (size_t gamma=0;gamma<int_dof;gamma++) {
			 	for (size_t gamma2=0;gamma2<int_dof;gamma2++) {
					Gtilde(i+gamma*n,j+gamma2*n) = setupGtildeAux(i,gamma,j,gamma2,G,dynVars,g,ether,aux);
				}
			}
		}
	}
				
}

void accOrbitalCorrelation(size_t gamma,size_t gamma2,
				Geometry const &geometry,DynVars const &dynVars, Parameters const &ether,Aux &aux)
{
        int n = geometry.volume();
        int int_dof =  int(ether.hilbertSize/n);
	
	size_t twoN  = 2*n;
	psimag::Matrix<MatType> G(twoN,twoN),Gtilde(twoN,twoN);	
	setupGreenFunction(G,ether,aux);
	setupGtilde(Gtilde,G,dynVars,geometry,ether,aux);
	
        // Do for all x's
	 for (size_t x=0;x<n;x++) {     // Sum over all i's (or w's)
                MatType tmp=0;
                for (size_t w=0;w<n;w++) {
                        //for (size_t gamma=0;gamma<int_dof;gamma++) {
                                //for (size_t gamma2=0;gamma2<int_dof;gamma2++) {
                                        size_t v=geometry.add(x,w);
                                        if (x==0 && gamma==gamma2) {
                                                tmp += (1.0-Gtilde(w+gamma*n,w+gamma*n));
                                        } else {
						tmp +=(1.0-Gtilde(w+gamma*n,w+gamma*n))*
						(1.0-Gtilde(v+gamma2*n,v+gamma2*n));
                                                tmp -= Gtilde(v+gamma2*n,w+gamma*n)*
						Gtilde(w+gamma*n,v+gamma2*n);
                                        }
                                //}
                        //}
                }
                aux.oco[x+gamma*n+gamma2*n*int_dof] += real(tmp)/n;
		// aux.cco[x] -= density*density; <-- ???
        }
	
}

//not correct, need "dressed" operators
void accOrbitalCorrelation2(Geometry const &geometry,DynVars const &dynVars, Parameters const &ether,Aux &aux)
{
	int i,j,k;
	int n = geometry.volume();
	MatType temp;
			
	for(i=0; i<n; i++) {
		for (j=0;j<n;j++) {
			k = geometry.add(i,j);
			temp += //calcOrbitalT(i,0,ether,aux)*calcOrbitalT(i,0,ether,aux)
			      //+ calcOrbitalT(i,1,ether,aux)*calcOrbitalT(i,1,ether,aux)
			      + calcOrbitalT(j,2,ether,aux)*calcOrbitalT(k,2,ether,aux);
		}
		aux.oco[i] += real(temp)/n;
	}
}

//! entry point for nanocluster (correlations)
template<typename FieldType>
void calcLocalk(psimag::Matrix<std::complex<FieldType> >& sq,const std::vector<size_t>& q,
		Geometry const &geometry,DynVars const &dynVars, Parameters const &ether)
{
	std::vector<FieldType> cd;
	std::vector<std::vector<size_t> > d;
	
	for (size_t plaquetteIndex=0;plaquetteIndex<ether.linSize;plaquetteIndex++) {
		calcCdAndD(plaquetteIndex,cd,d);
		calcSq(sq,q,cd,d,ether.kmesh,plaquetteIndex);
	}
}

template<typename FieldType>
FieldType calcCorrelation(size_t i,size_t j,bool doSpins,const DynVars& dynVars)
{
	if (!doSpins) 
		return dynVars.phonons[i][0]*dynVars.phonons[j][0]+
				dynVars.phonons[i][1]*dynVars.phonons[j][1]-
			(square(dynVars.phonons[i][0])+square(dynVars.phonons[i][1]));
	return cos(dynVars.theta[i])*cos(dynVars.theta[j])+
				sin(dynVars.theta[i])*sin(dynVars.theta[j])*
			cos(dynVars.phi[i]-dynVars.phi[j]);
	
}


template<typename DistanceType,typename FieldType>
void calcCdAndD(size_t plaquetteIndex,std::vector<FieldType>& cd,
		std::vector<DistanceType>& d,Geometry const &geometry,DynVars const &dynVars, Parameters const &ether)
{
	size_t n = ether.linSize;
	
	for (size_t i=0;i<n;i++) {
		for (size_t j=0;j<n;j++) {
			if (!geometry.isInPlaquette(plaquetteIndex,i)) continue;
			if (!geometry.isInPlaquette(plaquetteIndex,j)) continue;
 			DistanceType thisD;
			geometry.plaquetteDistance(thisD,i,j); // thisD = calcDistance(i,j);
			int x = isInVector(d,thisD);
			if (x<0) {
				d.push_back(thisD); // d[x] = thisD , 
				x = d.size()-1;	    // C(thisD) =  
				cd[x] = calcCorrelation<FieldType>(i,j,false,dynVars);
			} else  {
				cd[x] += calcCorrelation<FieldType>(i,j,false,dynVars);
			}
		}
	}
}


// q contains the indices for the k values we want to compute
template<typename FieldType,typename GeometryType>
void calcSq(psimag::Matrix<std::complex<FieldType> >& sq,,const std::vector<size_t>& q,const std::vector<FieldType>& cds,
	    	Kmesh& kmesh,size_t plaquetteIndex,const GeometryType& geometry)
{
	size_t nOfKs = q.size();
	for (size_t i=0;i<nOfKs;i++) { // loop over ks
		sq(i,plaquetteIndex) = std::complex<FieldType>(0,0);
		std::vector<size_t> tmp;
		kmesh.calcKVector(tmp,i); // put the i-th k-vector into tmp
		for (size_t j=0;j<cds.size();j++) { // loop over distances
			std::vector<size_t> dvector;
			geometry.plaquetteCalcD(j,dvector);
			
			FieldType factor = 2. * M_PI * scalarProduct(tmp,dvector)/kmesh.length();
			std::complex<FieldType> incr = std::complex<FieldType>(cos(factor),sin(factor)); // = exp(ik)
			sq(i,plaquetteIndex) += cds[j]*incr;
		}
		// FIXME: Is this normalization needed?
		//sqReal[i] /= $GlobalLc[0]*$GlobalLc[1];
		//$sqImag[$i] /= $GlobalLc[0]*$GlobalLc[1];
	}
}

