/*******************************************************************
 **         THIS FILE IS PART OF THE SPF (SPIN-PHONON-FERMION)    ** 
 **                          COMPUTER PROGRAM                     **
 **                                                               **
 ********************** SPF VERSION 6.2. ***************************
 **                                                               **
 **       For more info see: http://mri-fre.ornl.gov/spf          **
 ******************************************************************/

#include "common.h"
	
void auxCreateJmatrix(vector<MatType> &jmatrix, DynVars const &dynVars,Parameters const &ether,int site)
{
	int i,j;
	
        jmatrix[0]=0.5*cos(dynVars.theta[site]);

        jmatrix[1]=0.0;

        jmatrix[2]=MatType(0.5*sin(dynVars.theta[site])*cos(dynVars.phi[site])/sqrt(3.0),
                -0.5*sin(dynVars.theta[site])*sin(dynVars.phi[site])/sqrt(3.0));     
   
        jmatrix[3]=0.0;

        jmatrix[4]=MatType(sin(dynVars.theta[site])*cos(dynVars.phi[site])/sqrt(6.0),
                -sin(dynVars.theta[site])*sin(dynVars.phi[site])/sqrt(6.0));  

        jmatrix[5]=0.0;
        jmatrix[6]=0.0;

        jmatrix[7]=-cos(dynVars.theta[site])/6.0;

        jmatrix[8]=MatType(sin(dynVars.theta[site])*cos(dynVars.phi[site])/3.0,
                 sin(dynVars.theta[site])*sin(dynVars.phi[site])/3.0); 
     
        jmatrix[9]=MatType(0.5*sin(dynVars.theta[site])*cos(dynVars.phi[site])/sqrt(3.0),
                 -0.5*sin(dynVars.theta[site])*sin(dynVars.phi[site])/sqrt(3.0));

       jmatrix[10]=MatType(-sin(dynVars.theta[site])*cos(dynVars.phi[site])/(3.0*sqrt(2.0)),
                 -sin(dynVars.theta[site])*sin(dynVars.phi[site])/(3.0*sqrt(2.0)));
    
        jmatrix[11]=-2.0*cos(dynVars.theta[site])/sqrt(3.0);

        jmatrix[12]= MatType(0.5*sin(dynVars.theta[site])*cos(dynVars.phi[site])/sqrt(3.0),
                 0.5*sin(dynVars.theta[site])*sin(dynVars.phi[site])/sqrt(3.0)); 

         jmatrix[13]=MatType(sin(dynVars.theta[site])*cos(dynVars.phi[site])/3.0,
                 -sin(dynVars.theta[site])*sin(dynVars.phi[site])/3.0);

         jmatrix[14]= -jmatrix[7];
         jmatrix[15]=0.0;
         jmatrix[16]= jmatrix[11];
 
         jmatrix[17]=MatType(sin(dynVars.theta[site])*cos(dynVars.phi[site])/(3.0*sqrt(2.0)),
                 -sin(dynVars.theta[site])*sin(dynVars.phi[site])/(3.0*sqrt(2.0)));
    
          jmatrix[18]=0.0;

          jmatrix[19]=MatType(0.5*sin(dynVars.theta[site])*cos(dynVars.phi[site])/sqrt(3.0),
                 0.5*sin(dynVars.theta[site])*sin(dynVars.phi[site])/sqrt(3.0));   
  
         jmatrix[20]= 0.0;
         jmatrix[21]= -jmatrix[0];
         jmatrix[22]= 0.0;

         jmatrix[23]=MatType(-sin(dynVars.theta[site])*cos(dynVars.phi[site])/sqrt(6.0),
                -sin(dynVars.theta[site])*sin(dynVars.phi[site])/sqrt(6.0)); 

         jmatrix[24]=  MatType(sin(dynVars.theta[site])*cos(dynVars.phi[site])/sqrt(6.0),
                sin(dynVars.theta[site])*sin(dynVars.phi[site])/sqrt(6.0));  

         jmatrix[25]=MatType(-sin(dynVars.theta[site])*cos(dynVars.phi[site])/(3.0*sqrt(2.0)),
                 sin(dynVars.theta[site])*sin(dynVars.phi[site])/(3.0*sqrt(2.0)));

        jmatrix[26]= jmatrix[11];
        jmatrix[27]= 0.0;
        jmatrix[28]= jmatrix[7];

        jmatrix[29]= MatType(-sin(dynVars.theta[site])*cos(dynVars.phi[site])/6.0,
                 sin(dynVars.theta[site])*sin(dynVars.phi[site])/6.0);

         jmatrix[30]= 0.0;
         jmatrix[31]= jmatrix[11];

        jmatrix[32]=MatType(sin(dynVars.theta[site])*cos(dynVars.phi[site])/(3.0*sqrt(2.0)),
                 sin(dynVars.theta[site])*sin(dynVars.phi[site])/(3.0*sqrt(2.0)));

      jmatrix[33]=MatType(-sin(dynVars.theta[site])*cos(dynVars.phi[site])/sqrt(6.0),
                sin(dynVars.theta[site])*sin(dynVars.phi[site])/sqrt(6.0)); 

      jmatrix[34]= MatType(-sin(dynVars.theta[site])*cos(dynVars.phi[site])/6.0,
                 -sin(dynVars.theta[site])*sin(dynVars.phi[site])/6.0);
 
      jmatrix[35]= jmatrix[14];

	for (i=0;i<6;i++) {
		for (j=i+1;j<6;j++) {
			jmatrix[i+j*6]=conj(jmatrix[j+i*6]);
		}
	}
	
	for (i=0;i<36;i++) {
		jmatrix[i] *= ether.J[0];
		//cerr<<"i "<<jmatrix[i]<<endl;
	}

}


void kTpemHamiltonian (Geometry const &geometry, DynVars const &dynVars,
		 tpem_sparse *matrix,Parameters const &ether,Aux &aux,int type)
{
	int	row, col, volume, i = 0, p;
	double	a=aux.varTpem_a, b=aux.varTpem_b, tmp;
	int	j,k,dir;
	volume = ether.linSize;
	int level,gamma1,gamma2;
	vector<MatType> jmatrix(36,0);
	
	row=0;
	for (gamma1=0;gamma1<6;gamma1++) {
			for (p = 0; p < volume; p++, row++) {
			
			auxCreateJmatrix(jmatrix,dynVars,ether,p);
				
			matrix->rowptr[row] = i;
			matrix->colind[i] = row;		
			tmp =(real(jmatrix[gamma1+6*gamma1]))*double(ether.modulus[p]) + ether.potential[p] ;
			
			matrix->values[i] = (tmp-b)/a; /* -b because it's diagonal */
			i++;
			
				
			for (j = 0; j <  geometry.z(p); j++) {	
				k = geometry.neighbor(p,j);
				dir = geometry.scalarDirection(p,k,j);
  //                             for ($gamma1=0;$gamma1<ether.numberOfOrbitals;$gamma1++) {
				for (gamma2=0;gamma2<ether.numberOfOrbitals;gamma2++) {
					matrix->colind[i]=k+gamma2*volume;
					matrix->values[i]=ether.bandHoppings[gamma1+gamma2*6+36*dir]/a;
					i++;
				}
//			}
			
                     }
			
			for (gamma2=0;gamma2<6;gamma2++) {
				if (gamma1==gamma2) continue;
				matrix->colind[i]=p + gamma2*volume;
				matrix->values[i]=jmatrix[gamma1+6*gamma2]*double(ether.modulus[p])/a;
				i++;
			}
			
		}
	}
						
				
			
	if (i!=ether.nonzero) {
		cerr<<"i="<<i<<" but nonzero="<<ether.nonzero<<endl;
		cerr<<"At this point: "<<__FILE__<<" "<<__LINE__<<endl;
		exit(1);
	}
	
		
}

void setSupport(vector<unsigned int> &support,unsigned int i,Geometry const &geometry)
{
	int volume=geometry.volume();
	int shift;
	
	support.push_back(i);
	support.push_back(i+volume);
	support.push_back(i+2*volume);
	support.push_back(i+3*volume);
       	
	
}

void setHilbertParams(Parameters &ether, Aux &aux,Geometry const &geometry)
{
	int n=geometry.volume(), d=geometry.dim();
	int i;
	
	ether.typeofmodel="MODEL_KONDO_DMS_FCC";
	ether.hilbertSize=6*n;
	ether.nonzero= square(ether.numberOfOrbitals)*(1+geometry.z(0)) * n;
	
	if (!(geometry.latticeName()=="fcc")) {
		cerr<<"I was expecting fcc geometry but instead you input "<<geometry.latticeName()<<endl;
		cerr<<"At this point: "<<__FILE__<<" "<<__LINE__<<endl;
		exit(1);
	}
	
	cerr<<ether.bandHoppings.size()<<endl;
	vectorPrint(ether.bandHoppings,"hoppings=",cerr);
//	ether.energy1= -10-maxElement(ether.J);
        ether.energy1= -18.5;
	ether.energy2= 20.5;
	aux.varTpem_a = 0.5*(ether.energy2-ether.energy1);
	aux.varTpem_b = 0.5*(ether.energy2+ether.energy1);
	ether.classFieldList.push_back(0);
	
	
}
