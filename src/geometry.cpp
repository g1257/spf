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




#include "geometry.h"

/* SEE int main() at the end for an example of use */

/* Functions that should be exported 
  1. g_geometry(vector<siteIndex> &nn,vector<siteIndex> &nidx,int l,string const &type);
  type is: cubic, square, fcc, bcc, triangular  (hcp not ready yet)
  
  2. g_geometry(vector<siteIndex> &nn,vector<siteIndex> &nidx,vector<siteIndex> const &lvector,string const &type);
  type is: prism or rectangle
*/  
  
using namespace std;

/*
void split(std::string const &s,vector<siteIndex> &v)
{
	siteIndex i,l;
	std::string buffer="";
	
	l=s.size();
	for (i=0;i<l;i++) {
		if (s[i]==',' || s[i]=='x') {
			v.push_back(atoi(buffer.c_str()));
			buffer="";
		} else {
			buffer=buffer + s[i];
		}
	}
	v.push_back(atoi(buffer.c_str()));
	for (i=0;i<v.size();i++) {
		//cerr<<"split: Got v["<<i<<"]="<<v[i]<<endl;
	}
}
*/

bool Geometry::isCubicType(string const &s) const
{
	if (s=="square" || s=="cubic" ||
		s=="1d") return true;
	return false;
}

bool Geometry::isCubicType() const 
{
	return isCubicType(latticeType);
}


bool Geometry::g_pbc(int &x, int L) const
{
	bool r=false;
	if (x<0) r=true; 
	if (x>=L) r=true; 
	while(x<0) x+=L;
	while(x>=L) x-=L;
	return r;
}

int Geometry::g_index(int &x,int &y,int &z,int L)
{
	g_pbc(x,L);
	g_pbc(y,L);
	g_pbc(z,L);
	return x+y*L+z*L*L;
}

int Geometry::g_index(int &x,int &y,int &z,vector<siteIndex> const &lvector)
{
	g_pbc(x,lvector[0]);
	g_pbc(y,lvector[1]);
	if (lvector.size()>=3) g_pbc(z,lvector[2]);	
	
	return x+y*lvector[0]+z*lvector[0]*lvector[1];
}

/*
template<class T>
void Geometry::index2Coor(vector<T> &v,siteIndex i,std::string const &lt) const
{
	T x,y,z,l;
	l=Length[0];
	v.clear();
	if (lt=="1d") {
		v.push_back(i);
	} else if (lt=="square") {
		x=i % l;
		v.push_back(x);
		y=int(i/l);
		v.push_back(y);
	} else if (lt=="cubic") {
		z=int(i/(l*l));
		x=i % (l*l);
		y=int(x/l);
		x = x % l;
		v.push_back(x); v.push_back(y); v.push_back(z);
	} else {
		cerr<<"Geometry::index2Coor not implemented for latticeType="<<latticeType<<endl;
		cerr<<"At this point "<<__FILE__<<" "<<__LINE__<<endl;
		exit(1);
	}
} */



void Geometry::index2Coor(vector<double> &v,siteIndex i,std::string const &lt) const
{
	int x,y,z,l;
	l=Length[0];
	v.clear();
	if (lt=="1d") {
		v.push_back(i);
	} else if (lt=="square") {
		x=i % l;
		v.push_back(x);
		y=int(i/l);
		v.push_back(y);
	} else if (lt=="cubic") {
		z=int(i/(l*l));
		x=i % (l*l);
		y=int(x/l);
		x = x % l;
		v.push_back(x); v.push_back(y); v.push_back(z);
	} else {
		cerr<<"Geometry::index2Coor not implemented for latticeType="<<latticeType<<endl;
		cerr<<"At this point "<<__FILE__<<" "<<__LINE__<<endl;
		exit(1);
	}
} 

void Geometry::index2Coor(vector<int> &v,siteIndex i,std::string const &lt) const
{
	int x,y,z,l,j;
	l=Length[0];
	v.clear();
	if (lt=="1d") {
		v.push_back(i);
	} else if (lt=="square") {
		x=i % l;
		v.push_back(x);
		y=int(i/l);
		v.push_back(y);
	} else if (lt=="cubic") {
		z=int(i/(l*l));
		x=i % (l*l);
		y=int(x/l);
		x = x % l;
		v.push_back(x); v.push_back(y); v.push_back(z);
	} else if (lt=="prism") {
		j=(i-x)/Length[0];
		y = j % Length[1];
		z = int(j/Length[1]);
		x = i - y*Length[0] - z*Length[0] * Length[1];
		v.push_back(x); v.push_back(y); v.push_back(z);
	} else if (lt=="rectangle") {
		x = i % Length[0];
		y = int(i/Length[0]);
		v.push_back(x); v.push_back(y);
	} else	{
		cerr<<"Geometry::index2Coor not implemented for latticeType="<<latticeType<<endl;
		cerr<<"At this point "<<__FILE__<<" "<<__LINE__<<endl;
		exit(1);
	}
} 

/*template<class T>
void vectorPrint(vector<T> const &v,ostream &s)
{
	for (int i=0;i<v.size();i++) {
		s<<i<<" "<<v[i]<<endl;
	}
}*/

/*Sites are numbered in the following way: Bcc is considered as cubic plus a basis.
 * Sites of basis (0,0,0) are numbered first (first l*l*l) and those
 * with basis (a/2,a/2,a/2) are numbred next (last l*l*l).
 */
void Geometry::g_geometryBcc(vector<siteIndex> &nn,vector<siteIndex> &nidx,int l)
{
	int volume=l*l*l,x,y,z,xx,yy,zz,b,nbasis=2,shift;
	nn.insert(nn.begin(),volume*nbasis,8);
	for (b=0;b<nbasis;b++) {
		shift=2*b-1;
		for (z=0;z<l;z++) { for (y=0;y<l;y++) { for (x=0;x<l;x++) {
			xx=x; yy=y; zz=z;
			nidx.push_back(g_index(xx,yy,zz,l)+volume*(1-b));
			xx=x+shift;
			nidx.push_back(g_index(xx,yy,zz,l)+volume*(1-b));
			xx=x; yy=y+shift; zz=z;
			nidx.push_back(g_index(xx,yy,zz,l)+volume*(1-b));
			xx=x;yy=y;zz=z+shift;
			nidx.push_back(g_index(xx,yy,zz,l)+volume*(1-b));
			xx=x+shift;yy=y+shift;zz=z;
			nidx.push_back(g_index(xx,yy,zz,l)+volume*(1-b));
			xx=x+shift;yy=y;zz=z+shift;
			nidx.push_back(g_index(xx,yy,zz,l)+volume*(1-b));
			xx=x;yy=y+shift;zz=z+shift;
			nidx.push_back(g_index(xx,yy,zz,l)+volume*(1-b));
			xx=x+shift;yy=y+shift;zz=z+shift;			
			nidx.push_back(g_index(xx,yy,zz,l)+volume*(1-b));
		}}}
	}
}

/* See geometry_honeycomb.fig for numbering */
void Geometry::g_geometryHoneycomb(vector<siteIndex> &nn,vector<siteIndex> &nidx,vector<int> const &lvector)
{
	int volume=lvector[0]*lvector[1],x,y,tmp,b,nbasis=2,shift,type;
	nn.insert(nn.begin(),volume*nbasis,3);
	
	for (type=0;type<2;type++) {
		for (y=0;y<lvector[1];y++) {
			for (x=0;x<lvector[0];x++) {
				
				nidx.push_back(x+y*lvector[0]+(1-type)*volume);
				border.push_back(-1);
				
				tmp = x-2*type+1;
				b=g_pbc(tmp,lvector[0]);
				if (b) border.push_back(0);
				else border.push_back(-1);
				nidx.push_back(tmp+y*lvector[0]+(1-type)*volume);
				
				if (y==0 || y%2==0) tmp = y+2*type-1;
				else tmp = y-2*type+1;
				b=g_pbc(tmp,lvector[1]);
				if (b) border.push_back(1);
				else border.push_back(-1);
				nidx.push_back(x+tmp*lvector[0]+type*volume);
			}
		}
	}
}
	
				
			
void Geometry::g_add(int &xx,int &yy, int &zz,vector<int> const &v) const
{
	xx+=v[0];
	yy+=v[1];
	zz+=v[2];
}

void Geometry::g_add(vector<double> &dest,vector<double> const &src) const
{
	int i;
	for (i=0;i<src.size();i++) dest[i] += src[i];
	
}


double Geometry::g_norm(vector<double> const &r)
{
	int i;
	double s=0;
	for (i=0;i<r.size();i++) s+= r[i]*r[i];
	return s;
}


double Geometry::g_auxZb_distance(vector<double> const &r1,vector<double> const &r2,g_aux_Zb_data const &data)
{
	int l=data.length;
	int i;
	double minDistance=data.minDistance;
	vector<double> r(r1.size());
	
	for (i=0;i<r1.size();i++) {
		r[i]=r1[i]-r2[i];
		while(r[i]<0.25) r[i]+=l;
		while(r[i]>=l-0.25) r[i]-=l;
	}
	return g_norm(r);

}

void Geometry::g_auxZb_initBasis(vector<vector<double> > &basisVector) const
{
	vector<double> tmpVector;
	
	
	tmpVector.push_back(0); tmpVector.push_back(0); tmpVector.push_back(0);
	basisVector.push_back(tmpVector);
	
	tmpVector[0]=0; tmpVector[1]=0.5; tmpVector[2]=0.5;
	basisVector.push_back(tmpVector);
	
	tmpVector[0]=0.5; tmpVector[1]=0; tmpVector[2]=0.5;
	basisVector.push_back(tmpVector);
	
	tmpVector[0]=0.5; tmpVector[1]=0.5; tmpVector[2]=0;
	basisVector.push_back(tmpVector);
	
	tmpVector[0]=0.25; tmpVector[1]=0.25; tmpVector[2]=0.25;
	basisVector.push_back(tmpVector);
	
	tmpVector[0]=0.25; tmpVector[1]=0.75; tmpVector[2]=0.75;
	basisVector.push_back(tmpVector);
	
	tmpVector[0]=0.75; tmpVector[1]=0.25; tmpVector[2]=0.75;
	basisVector.push_back(tmpVector);
	
	tmpVector[0]=0.75; tmpVector[1]=0.75; tmpVector[2]=0.25;
	basisVector.push_back(tmpVector);
	
	
	
}




// order of the 8 basis vectors is
// (0,0,0) (0,a/2,a/2), (a/2,0,a/2), (a/2,a/2,0), and these + (a/4,a/4,a/4)
void Geometry::g_auxZb_Index2Coor(vector<double> &r,int site,int volume) const
{
	// coordinates are (x,y,z,b) where b is the basis
	static int firstcall=0;
	static vector<vector<double> > basisVector;
	
	if (firstcall==0) {
		firstcall=1;
		g_auxZb_initBasis(basisVector);
	}
	
	int b = int(site/volume);
	siteIndex i = site % volume;
	std::string s="cubic";
	index2Coor(r,i,s);
	g_add(r,basisVector[b]);
}
	


double Geometry::g_auxZb_dist(siteIndex site1, siteIndex site2,g_aux_Zb_data const &data)
{
	vector<double> r1,r2;
	
	g_auxZb_Index2Coor(r1,site1,data.volume); // ordering dependence
	g_auxZb_Index2Coor(r2,site2,data.volume); // ordering dependence 
	return g_auxZb_distance(r1,r2,data);
}

/* Could be too slow */
void Geometry::g_auxZb_findNn(vector<siteIndex> &v,int site,g_aux_Zb_data const &data)
{
	int i;
	double dd;
	v.clear();
	for (i=0;i<data.volume*data.nbasis;i++) {
		dd=g_auxZb_dist(i,site,data);
		if (fabs(dd-data.minDistance)<1e-6) {
			v.push_back(i);
		}
	}
	if (v.size()!=Nn[site]) {
		cerr<<"Site "<<site<<" has "<<v.size()<<" neighbours but it should be "<<Nn[site]<<" instead\n";
		exit(1);
	}
}

void Geometry::g_geometryZb(vector<siteIndex> &nn,vector<siteIndex> &nidx,int l)
{
	int i,j,vol=l*l*l,nbasis=8;
	vector<siteIndex> tmpVec;
	nn.insert(nn.begin(),vol*nbasis,4);
	g_aux_Zb_data data;
	data.minDistance=3.0/16.0;
	data.length=l;
	data.nbasis=nbasis;
	data.volume=vol;
	
	for (i=0;i<vol*nbasis;i++) {
		g_auxZb_findNn(tmpVec,i,data); // find sites that are a distance min from i and store them on tmpVec
		for (j=0;j<tmpVec.size();j++) nidx.push_back(tmpVec[j]);
	}
}

void Geometry::g_geometryZbDirection(vector<double> &r,siteIndex site1, siteIndex site2) const
{
	vector<double> r1,r2;
	int i;
	int l=Length[0];
	
	
	g_auxZb_Index2Coor(r1,site1,l*l*l);
	// cerr<<"HEREEEEEEEEEEEEE "<<r1.size()<<" "<<endl;
	g_auxZb_Index2Coor(r2,site2,l*l*l);
	
	r.clear();
	//cerr<<"Direction site "<<site1<<" from site "<<site2<<endl;
	
	for (i=0;i<r1.size();i++) {
		//cerr<<r1[i]<<" "<<r2[i]<<endl;
		r.push_back(r1[i]-r2[i]);
		//cerr<<r[i]<<" ";
		if (r[i]< -0.25) r[i]+=l;
		if (r[i]>=l-0.25) r[i]-=l;
	}
	//cerr<<endl;
}

/*************** FCC *****************************/
void Geometry::g_auxFcc_initBasis(vector<vector<double> > &basisVector) const
{
	vector<double> tmpVector;
	
	
	tmpVector.push_back(0); tmpVector.push_back(0); tmpVector.push_back(0);
	basisVector.push_back(tmpVector);
	
	tmpVector[0]=0; tmpVector[1]=0.5; tmpVector[2]=0.5;
	basisVector.push_back(tmpVector);
	
	tmpVector[0]=0.5; tmpVector[1]=0; tmpVector[2]=0.5;
	basisVector.push_back(tmpVector);
	
	tmpVector[0]=0.5; tmpVector[1]=0.5; tmpVector[2]=0;
	basisVector.push_back(tmpVector);

	
}

int getBasisIndex(vector<vector<double> > const &basisVectors,vector<double> const &v)
{
	int i;
	for (i=0;i<basisVectors.size();i++) {
		if (v==basisVectors[i]) return i;
	}
	return -1; 
}

// order of the 4 basis vectors is
// (0,0,0) (0,a/2,a/2), (a/2,0,a/2), (a/2,a/2,0)
void Geometry::g_auxFcc_Index2Coor(vector<double> &r,int site,int vol) const
{
	// coordinates are (x,y,z,b) where b is the basis
	static int firstcall=0;
	static vector<vector<double> > basisVector;
	
	if (firstcall==0) {
		firstcall=1;
		g_auxFcc_initBasis(basisVector);
	}
	
	
	int b = int(site/vol);
	siteIndex i = site % vol;
	std::string s="cubic";
	index2Coor(r,i,s);
	g_add(r,basisVector[b]);
	
}

double Geometry::g_auxFcc_distance(vector<double> const &r1,vector<double> const &r2,g_aux_Zb_data const &data)
{
	int l=data.length;
	int i;
	double minDistance=data.minDistance;
	vector<double> r(r1.size());
	
	for (i=0;i<r1.size();i++) {
		r[i]=r1[i]-r2[i];
		while(r[i]<0.5) r[i]+=l;
		while(r[i]>=l-0.5) r[i]-=l;
	}
	return g_norm(r);

}

double Geometry::g_auxFcc_dist(siteIndex site1, siteIndex site2,g_aux_Zb_data const &data)
{
	vector<double> r1,r2;
	int l = Length[0];
	
	g_auxFcc_Index2Coor(r1,site1,l*l*l);
	g_auxFcc_Index2Coor(r2,site2,l*l*l);
	return g_auxFcc_distance(r1,r2,data);
}

void Geometry::g_auxFcc_findNn(vector<siteIndex> &v,int site,g_aux_Zb_data const &data)
{
	int i;
	double dd;
	v.clear();
	for (i=0;i<data.volume*data.nbasis;i++) {
		dd=g_auxFcc_dist(i,site,data);
		if (fabs(dd-data.minDistance)<1e-6) {
			v.push_back(i);
		}
	}
	if (v.size()!=Nn[site]) {
		cerr<<"(FCC) Site "<<site<<" has "<<v.size()<<" neighbours but it should be "<<Nn[site]<<" instead\n";
		exit(1);
	}
}
	
void Geometry::g_geometryFcc(vector<siteIndex> &nn,vector<siteIndex> &nidx,int l)
{
	int i,j,vol=l*l*l,nbasis=4;
	vector<siteIndex> tmpVec;
	nn.insert(nn.begin(),vol*nbasis,12);
	g_aux_Zb_data data;
	data.minDistance=0.5;
	data.length=l;
	data.nbasis=nbasis;
	data.volume=vol;
	
	for (i=0;i<vol*nbasis;i++) {
		g_auxFcc_findNn(tmpVec,i,data); // find sites that are a distance min from i and store them on tmpVec
		for (j=0;j<tmpVec.size();j++) nidx.push_back(tmpVec[j]);
	}
}

void Geometry::g_geometryFccDirection(vector<double> &r,siteIndex site1, siteIndex site2) const
{
	vector<double> r1,r2;
	int i;
	int l=Length[0];
		
	g_auxFcc_Index2Coor(r1,site1,l*l*l);
	g_auxFcc_Index2Coor(r2,site2,l*l*l);
	
	r.clear();
	for (i=0;i<r1.size();i++) {
		r.push_back(r1[i]-r2[i]);
		while(r[i]<0.5) r[i]+=l;
		while(r[i]>=l-0.5) r[i]-=l;
		
	}
}

void Geometry::direction(vector<double> &r,siteIndex site1,siteIndex site2) const
{
	if (latticeType=="fcc") {
		g_geometryFccDirection(r,site1,site2);
	} else if (latticeType=="zincblende") {
		g_geometryZbDirection(r,site1,site2);
	} else {
		cerr<<"DIRECTION NOT IMPLEMENTED\n";
		exit(1);
	}
}
/*
int Geometry::scalarDirectionOld(siteIndex site1,siteIndex site2,int number) const
{
	vector<double> r(3);
//	int i;
	if (number>=0 && isCubicType()) return int(number/2);
	
	direction(r,site1,site2);

          if (r[0]<0) {
                r[0] = -r[0];
                r[1] = -r[1];
                r[2] = -r[2];
        }
//         if (r[0]==0.25 && r[1]==0.25 && r[2]==0.25) return 0;
//         if (r[0]==0.25 && r[1]==-0.25 && r[2]==-0.25) return 1;
//         if (r[0]==0.25 && r[1]==-0.25 && r[2]==0.25) return 2;
//         if (r[0]==0.25 && r[1]==0.25 && r[2]==-0.25) return 3;
          vector<int> r2(3); 
          int i;    
          for(i=0;i<3;i++) r2[i]=int(r[i]*4);      

	  if (r2[0]==1 && r2[1]==1 && r2[2]==1) return 0;
          if (r2[0]==1 && r2[1]==-1 && r2[2]==-1) return 1;      
          if (r2[0]==1 && r2[1]==-1 && r2[2]==1) return 2;
          if (r2[0]==1 && r2[1]==1 && r2[2]==-1) return 3;
	 
//	cerr<<"SCALAR DIRECTION NOT IMPLEMENTED\n";
//	exit(1);
        std::cerr<<"Vectorial direction "<<r2[0]<<" "<<r2[1]<<" "<<r2[2];
        std::cerr<<" does not have a corresponding scalar direction.\n";
        std::cerr<<"At this point: "<<__FILE__<<" "<<__LINE__<<endl;

        exit(1);
        return -1; // unreachable, to avoid compiler warnings.
      
        
}
 */
 
 
void scDirAux(vector<double> const &r, int &tmpx, int &tmpy, int &tmpl)
{
	if (r[0]==0) {
		tmpl=0;
		tmpx=r[1]+0.5;
		tmpy=r[2]+0.5;
	} else if (r[1]==0) {
		tmpl=1;
		tmpx=r[0]+0.5;
		tmpy=r[2]+0.5;
	} else if (r[2]==0) {
		tmpl=2;
		tmpx=r[0]+0.5;
		tmpy=r[1]+0.5;
	} else {
		cerr<<"Error: vector r has no entry 0: "<<r[0]<<" "<<r[1]<<" "<<r[2]<<endl;
	}
}	
		
int Geometry::scalarDirection(siteIndex site1,siteIndex site2,int number) const
{
	vector<double> r(3);
//	int i;
	if (number>=0 && isCubicType()) return int(number/2);
	int tmpx,tmpy,tmpl;
	
	direction(r,site1,site2);


          vector<int> r2(3); 
          int i;    
if (latticeType=="zincblende") {
          for(i=0;i<3;i++) r2[i]=int(r[i]*4);      

	  if (r2[0]==1 && r2[1]==1 && r2[2]==1) return 0;
          if (r2[0]==1 && r2[1]==-1 && r2[2]==-1) return 1;      
          if (r2[0]== -1 && r2[1]==1 && r2[2]== -1) return 2;
          if (r2[0]== -1 && r2[1]== -1 && r2[2]== 1) return 3;
	 
	 for(i=0;i<3;i++) r2[i] = -r2[i];
	 
	 if (r2[0]==1 && r2[1]==1 && r2[2]==1) return 4;
          if (r2[0]==1 && r2[1]==-1 && r2[2]==-1) return 5;      
          if (r2[0]== -1 && r2[1]==1 && r2[2]== -1) return 6;
          if (r2[0]== -1 && r2[1]== -1 && r2[2]== 1) return 7;
} else if (latticeType=="fcc") {
	
	
	scDirAux(r,tmpx,tmpy,tmpl);
	return tmpx+2*tmpy+4*tmpl;	
}		  
//	cerr<<"SCALAR DIRECTION NOT IMPLEMENTED\n";
//	exit(1);
        std::cerr<<"Vectorial direction "<<r2[0]<<" "<<r2[1]<<" "<<r2[2];
        std::cerr<<" does not have a corresponding scalar direction.\n";
        std::cerr<<"At this point: "<<__FILE__<<" "<<__LINE__<<endl;

        exit(1);
        return -1; // unreachable, to avoid compiler warnings.
      
        
} 

/*********************** TRIANGULAR ********************************/
void Geometry::g_geometryTriangular(vector<siteIndex> &nn,vector<siteIndex> &nidx,int l)
{
	int volume=l*l,x,y,xx,yy,zz,b,nbasis=2;
	nn.insert(nn.begin(),volume*nbasis,6);
	zz=0;
	for (b=0;b<nbasis;b++) {
		for (y=0;y<l;y++) { for (x=0;x<l;x++) {
			xx=x+1; yy=y;;
			nidx.push_back(g_index(xx,yy,zz,l)+volume*b);
			xx=x-1; yy=y;
			nidx.push_back(g_index(xx,yy,zz,l)+volume*b);
			xx=x; yy=y;
			nidx.push_back(g_index(xx,yy,zz,l)+volume*(1-b));
			xx=x+2*b-1; yy=y;
			nidx.push_back(g_index(xx,yy,zz,l)+volume*(1-b));
			xx=x; yy=y+1;
			nidx.push_back(g_index(xx,yy,zz,l)+volume*(1-b));
			xx=x-1; yy=y+1;
			nidx.push_back(g_index(xx,yy,zz,l)+volume*(1-b));
		}}
	}
}


void Geometry::g_geometryPrism(vector<siteIndex> &nn,vector<siteIndex> &nidx,vector<siteIndex> const &lvector)
{
	int volume=lvector[0]*lvector[1]*lvector[2],x,y,z,xx,yy,zz;
	nn.insert(nn.begin(),volume,6);
	for (z=0;z<lvector[2];z++) { for (y=0;y<lvector[1];y++) { for (x=0;x<lvector[0];x++) {		
		xx=x+1; yy=y; zz=z;
		if (g_pbc(xx,lvector[0])) {
			border.push_back(0);
		} else {
			border.push_back(-1);
		}
		nidx.push_back(g_index(xx,yy,zz,lvector));
		xx=x-1;
		if (g_pbc(xx,lvector[0])) {
			border.push_back(0);
		} else {
			border.push_back(-1);
		}
		nidx.push_back(g_index(xx,yy,zz,lvector));
		xx=x; yy=y+1;
		if (g_pbc(yy,lvector[1])) {
			border.push_back(1);
		} else {
			border.push_back(-1);
		}
		nidx.push_back(g_index(xx,yy,zz,lvector));
		yy=y-1;
		if (g_pbc(yy,lvector[1])) {
			border.push_back(1);
		} else {
			border.push_back(-1);
		}
		nidx.push_back(g_index(xx,yy,zz,lvector));
		yy=y; zz=z+1;
		if (g_pbc(zz,lvector[2])) {
			border.push_back(2);
		} else {
			border.push_back(-1);
		}
		nidx.push_back(g_index(xx,yy,zz,lvector));
		zz=z-1;
		if (g_pbc(zz,lvector[2])) {
			border.push_back(2);
		} else {
			border.push_back(-1);
		}
		nidx.push_back(g_index(xx,yy,zz,lvector));
	}}}
}

void Geometry::g_geometryCubic(vector<siteIndex> &nn,vector<siteIndex> &nidx,int l)
{
	vector<siteIndex> lvector(3,l);
	g_geometryPrism(nn,nidx,lvector);
}


void Geometry::g_geometryRectangle(vector<siteIndex> &nn,vector<siteIndex> &nidx,vector<siteIndex> const &lvector)
{
	int volume=lvector[0]*lvector[1],x,y,z,xx,yy,zz=0;
	nn.insert(nn.begin(),volume,4);
	for (y=0;y<lvector[1];y++) { for (x=0;x<lvector[0];x++) {		
		xx=x+1; yy=y;
		if (g_pbc(xx,lvector[0])) {
			border.push_back(0);
		} else {
			border.push_back(-1);
		}
		nidx.push_back(g_index(xx,yy,zz,lvector));
		xx=x-1;
		if (g_pbc(xx,lvector[0])) {
			border.push_back(0);
		} else {
			border.push_back(-1);
		}
		nidx.push_back(g_index(xx,yy,zz,lvector));
		xx=x; yy=y+1;
		if (g_pbc(yy,lvector[1])) {
			border.push_back(1);
		} else {
			border.push_back(-1);
		}
		nidx.push_back(g_index(xx,yy,zz,lvector));
		yy=y-1;
		if (g_pbc(yy,lvector[1])) {
			border.push_back(1);
		} else {
			border.push_back(-1);
		}
		nidx.push_back(g_index(xx,yy,zz,lvector));
	}}
}

void Geometry::g_geometryOneD(vector<siteIndex> &nn,vector<siteIndex> &nidx,int l)
{
	int x,xx;
	
	for (x=0;x<l;x++) {
		nn.push_back(2);
		xx=x+1;
		if (g_pbc(xx,l)) {
			border.push_back(0);
		} else {
			border.push_back(-1);
		}
		nidx.push_back(xx);
		xx=x-1;
		if (g_pbc(xx,l)) {
			border.push_back(0);
		} else {
			border.push_back(-1);
		}
		nidx.push_back(xx);	
	}
}

bool Geometry::allSidesEqual() const
{
	if (latticeType=="square" || latticeType=="cubic" || latticeType=="fcc" || 
		latticeType=="1d" || latticeType=="triangular" || latticeType=="zincblende") return true;
	return false;
}



void Geometry::g_geometryDistance2(vector<siteIndex> &nn2,vector<siteIndex> &nidx2,
	vector<siteIndex> const &lvector)
{
	siteIndex x,y,xx,yy,zz;
	
	if (!(latticeType=="square")) {
		cerr<<"Geometry::g_geometryDistance2 not implemented for latticeType="<<latticeType<<endl;
		exit(1);
	}
	zz=0;
	for (y=0;y<lvector[1];y++) { for (x=0;x<lvector[0];x++) {	
		nn2.push_back(4);	
		xx=x+1; yy=y+1;
		nidx2.push_back(g_index(xx,yy,zz,lvector));
		xx=x-1; yy=yy-1;
		nidx2.push_back(g_index(xx,yy,zz,lvector));
		xx=x-1; yy=y+1;
		nidx2.push_back(g_index(xx,yy,zz,lvector));
		xx=x+1; yy=y-1;
		nidx2.push_back(g_index(xx,yy,zz,lvector));
	}}
}



// Generalize: FIXME NONCUBIC NON-CUBIC (supports prism and rectangle)
siteIndex Geometry::coor2Index(vector<siteIndex> &v,string const &lt) const
{
	siteIndex pos,temp;
	siteIndex D = v.size();
	siteIndex i;
	
	if (!isCubicType(lt) && !(lt=="prism") && !(lt=="rectangle")) {
		cerr<<"Geometry::coor2Index not implemented for latticeType="<<latticeType<<endl;
		cerr<<"At this point "<<__FILE__<<" "<<__LINE__<<endl;
		exit(1);
	}
	g_pbc(v[0],Length[0]);
	pos=v[0];
	temp=1;

	for (i=1;i<D;i++) {
		g_pbc(v[i],Length[i]);
		temp*=Length[i-1];
		pos+=v[i]*temp;
	}	
	return pos;
}

// Generalize: FIXME NONCUBIC NON-CUBIC
double Geometry::scalarProd(siteIndex i,siteIndex j) const 
{
	double sum=0;
	vector<siteIndex> tmp,tmp2;
	
	if (!isCubicType()) {
		cerr<<"Geometry::scalarProd not implemented for latticeType="<<latticeType<<endl;
		cerr<<"At this point "<<__FILE__<<" "<<__LINE__<<endl;
		return 0.0;
		//exit(1);
	}
	
	index2Coor(tmp,i,latticeType);
	index2Coor(tmp2,j,latticeType);
	for (i=0;i<tmp.size();i++) sum+= tmp[i]*tmp2[i];
			
	return sum;
				
}

int Geometry::add_fcc_aux(vector<double> const &b1,vector<double> const &b2) const
{
	unsigned int i;
	int flaghalf=0;
	int flagzero=0;
	vector<double> b(b1.size());
	b= b1+b2;
	for (i=0;i<b.size();i++) {
		if (b[i]==0.5) flaghalf=1;
		if (b[i]==0) flagzero=1;
	}
	if (flaghalf==0 && flagzero==1) return  0;
	if (flaghalf==1 && flagzero==0) return 1;
	if (flaghalf==1 && flagzero==1) return 2;
	cerr<<"Geometry::add_fcc_aux: INTERNAL ERROR.\n";
	cerr<<"Exiting at this point "<<__FILE__<<" "<<__LINE__<<"\n";
	exit(1); //ERROR CODE
	return -1; //to satisfy the compiler	
}

// output rsum + bsum = input r1 + b1 +r2 +b2
void Geometry::add_fcc(vector<double> const &r1,vector<double> const &b1,vector<double> const &r2,vector<double> const &b2,
vector<double> &rsum,vector<double> &bsum) const 
{
	vector<double> unityvector;
	unityvector.push_back(1); unityvector.push_back(1); unityvector.push_back(1);
	int casetype=add_fcc_aux(b1,b2);
	
	switch (casetype) {
		case 0:  // b1+b2 == (a,a,0) or permutations
			rsum = r1+r2+b1+b2;
			vectorSet(bsum,0);
			break;
		case 1: // b1 + b2 == (a/2,a/2,a) or permutations
			bsum = unityvector-b1 -b2;
			rsum = r1 + r2 + 2*(b1+b2)-unityvector;
			break;
		case 2: // b1 + b2 == (a/2,a/2,0) or permutations
			rsum =r1+r2;
			bsum = b1+b2;
			break;
	}
		
}

int Geometry::add_fcc(siteIndex ind,siteIndex ind2) const
{
	// coordinates are (x,y,z,b) where b is the basis
	static int firstcall=0;
	static vector<vector<double> > basisVector;
	std::string s="cubic";
	int indr,ind2r,b,b2,vol=Length[0]*Length[0]*Length[0];
	vector<double> rvector(3),r2vector(3),rsumvector(3),bsumvector(3);
	vector<siteIndex> r(3);
	siteIndex i;
	
	if (firstcall==0) {
		firstcall=1;
		g_auxFcc_initBasis(basisVector);
	}
		
	// transform ind --> rvector, bvector and ind2 --> r2vector, b2vector
	b = int(ind/vol);
	indr = ind % vol;
	index2Coor(rvector,indr,s);
	b2 = int(ind2/vol);
	ind2r= ind2 % vol;
	index2Coor(r2vector,ind2r,s);
	//vectorPrint(rvector,"rvector");
	//vectorPrint(basisVector[b],"bvector");
	//vectorPrint(r2vector,"r2vector");
	//vectorPrint(basisVector[b2],"b2vector");
	
	// sum these two and obtain rsumvector bsumvector
	add_fcc(rvector,basisVector[b],r2vector,basisVector[b2],rsumvector,bsumvector);
	//vectorPrint(rsumvector,"rsumvector");
	//vectorPrint(bsumvector,"bsumvector");
	
	// transform rsumvector and bsumvector into an index that is returned
	b = getBasisIndex(basisVector,bsumvector);
	if (b<0) {
		cerr<<"Geometry::add (Geometry::add_fcc): INTERNAL ERROR.\n";
		cerr<<"Exiting at this point "<<__FILE__<<" "<<__LINE__<<"\n";
		exit(1); //ERROR CODE
	}
	for (i=0;i<rsumvector.size();i++) r[i] = rsumvector[i];
	indr = coor2Index(r,s);
	return (indr + b*vol);
}


// Generalize: FIXME NONCUBIC NON-CUBIC (added for prism and rectangle)
int Geometry::add(siteIndex ind,siteIndex ind2) const
{
	siteIndex i;
	vector<siteIndex> x,y;
	siteIndex l=Length[0];
	
	
	if (latticeType=="fcc") return add_fcc(ind,ind2);
	
	if (!isCubicType() && !(latticeType=="prism") && !(latticeType=="rectangle")) {
		if (errorPrinted==0) {
			cerr<<"Geometry::add not implemented for latticeType="<<latticeType<<endl;
			cerr<<"At this point "<<__FILE__<<" "<<__LINE__<<endl;
			errorPrinted=1;
		}
		return 0;
		//exit(1);
	}
	
	
	index2Coor(x,ind,latticeType);
	index2Coor(y,ind2,latticeType);
	//cerr<<"ind="<<ind<<" and ind2="<<ind2<<endl;
	for (i=0;i<x.size();i++) {
		//cerr<<"add: i="<<i<<" x="<<x[i]<<" y="<<y[i]<<endl;
		x[i] += y[i];
		g_pbc(x[i],Length[i]);
		//cerr<<"add: after g_pbc: x["<<i<<"]="<<x[i]<<endl;
	}
	i=coor2Index(x,latticeType);
	return i;
}

// Generalize: FIXME NONCUBIC NON-CUBIC
int Geometry::dim() const
{
	if (latticeType=="cubic" || latticeType=="prism" || latticeType=="fcc" || latticeType=="zincblende") return 3;
	if (latticeType=="square" || latticeType=="rectangle" || latticeType=="honeycomb") return 2;
	if (latticeType=="1d") return 1;
	cerr<<"Geometry::dim(): Not all lattice sides equal (not implemented)\n";
	exit(1);
	return 0; // to satisfy the compiler 
}

// Generalize: FIXME NONCUBIC NON-CUBIC
int Geometry::length() const
{
	
	if (!allSidesEqual()) {
		cerr<<"Geometry::length not implemented for latticeType="<<latticeType<<endl;
		exit(1);
	}
	return Length[0]; // to satisfy the compiler 
}

void Geometry::g_geometrySquare(vector<siteIndex> &nn,vector<siteIndex> &nidx,int l)
{
	vector<siteIndex> lvector(2,l);
	g_geometryRectangle(nn,nidx,lvector);
	g_geometryDistance2(Nn2,Nidx2,lvector);
}	

void Geometry::geometry(vector<siteIndex> &nn,vector<siteIndex> &nidx,int l,string const &type)
{
	siteIndex i;
	vector<int> lvector;
	
	if (type=="cubic") {
		g_geometryCubic(nn,nidx,l);
		for (i=0;i<3;i++) Length.push_back(l);
	} else if (type=="bcc") {
		g_geometryBcc(nn,nidx,l);
		for (i=0;i<3;i++) Length.push_back(l);
	} else if (type=="fcc") {
		for (i=0;i<3;i++) Length.push_back(l);
		g_geometryFcc(nn,nidx,l);
	} else if (type=="triangular") {
		g_geometryTriangular(nn,nidx,l);
		for (i=0;i<2;i++) Length.push_back(l);
	} else if (type=="square") {
		g_geometrySquare(nn,nidx,l);
		for (i=0;i<2;i++) Length.push_back(l);
	} else if (type=="1d") {
		g_geometryOneD(nn,nidx,l);
		for (i=0;i<1;i++) Length.push_back(l);
	} else if (type=="zincblende") {
		for (i=0;i<3;i++) Length.push_back(l);
		g_geometryZb(nn,nidx,l);
	} else {
		cerr<<"Unknown lattice type: *"<<type<<"*"<<endl;
		cerr<<"At this point "<<__FILE__<<" "<<__LINE__<<endl;
		exit(1);
	}
}

void Geometry::geometry(vector<siteIndex> &nn,vector<siteIndex> &nidx,vector<siteIndex> const &lvector,
string const &type)
{
	if (type=="prism") {
		g_geometryPrism(nn,nidx,lvector);
	} else if (type=="rectangle") {
		g_geometryRectangle(nn,nidx,lvector);
	} else if (type=="honeycomb") {
		g_geometryHoneycomb(nn,nidx,lvector);
	} else {
		cerr<<"Unknown lattice type: *"<<type<<"*"<<endl;
		cerr<<"At this point "<<__FILE__<<" "<<__LINE__<<endl;
		exit(1);
	}
}


void Geometry::nicePrint(ostream &s,int option) const
{
	int i,j,counter=0,tmp;
	vector<double> r(3);
	
	/* g_auxZb_Index2Coor(r,32,Length[0]*Length[0]*Length[0]);
	for (i=0;i<r.size();i++) cerr<<r[i]<<" ";
	cerr<<endl;
	exit(1); */ 
	
	for (i=0;i<Nn.size();i++) {
		for (j=0;j<Nn[i];j++) {
			s<<"The "<<j<<"-th neighbour of site "<<i<<" is "<<Nidx[counter];
			if (option) s<<" with border="<<borderId(i,j);
			if (latticeType == "zincblende" ||  latticeType =="fcc") {
				//r=pdirection(i,Nidx[counter]);
				direction(r,i,Nidx[counter]);
				s<<" direction is "<<r[0]<<" "<<r[1]<<" "<<r[2];
				tmp=scalarDirection(i,Nidx[counter],j);
				s<<" scalardir="<<tmp<<endl;
				
			}
			s<<endl;
			counter++;
		}
	}
}

void Geometry::init(string const &type,vector<siteIndex> const &Lvector)
{
	//while (Length.size()<3) Length.push_back(1);
	Length=Lvector;
	latticeType=type;
	//while (Length.size()<3) Length.push_back(1);
	geometry(Nn,Nidx,Length,latticeType);
	isInit=true;
}

void Geometry::init(string const &type,int l)
{
	
	latticeType=type;
	geometry(Nn,Nidx,l,latticeType);
	isInit=true;
}
/*!\brief Initializer that takes the lattice type and a comma separated list of lengths.
	
	* \param s1 can be "1d", "square", "cubic", "triangular" (these for have all sides equal)
	*          or "honeycomb", "prism" or "rectangular" (sides can be different)
	* \param s2 can be L for "1d", L,L for "square" or "triangular", L,L,L for "cubic",
	*        L1xL2 for "rectangular" or "honeycomb" or L1xL2xL3 for "prism"
	*/
void Geometry::init(string const &s1,string const &s2,size_t side,int verbose)
{
	vector<siteIndex> lvector;
	mysplit(s2,lvector,',');
	latticeType=s1;
	if (allSidesEqual()) init(s1,lvector[0]);
	else init(s1,lvector);
	if (verbose) nicePrint(cerr);
	plaquette_.init(this,side);
}
		
int Geometry::z(int index,int distance) const
{
	if (index<0 || index>=volume()) {
		cerr<<"Geometry::z: index="<<index<<" out of range\n";
		cerr<<"At this point: "<<__FILE__<<" "<<__LINE__<<endl;
		exit(1);
	}
	if (distance==1) return Nn[index];
	if (!(latticeType=="square") || distance>2) {
		cerr<<"Geometry::z: distance>2 or distance==2 && type!=square\n";
		cerr<<"             (not implemented)\n";
		cerr<<"At this point: "<<__FILE__<<" "<<__LINE__<<endl;
		exit(1);
	}
	return Nn2[index];
}

int Geometry::neighbor(int index,int n,int distance) const 
{
	if (index<0 || index>=volume()) {
		cerr<<"Geometry::z: index="<<index<<" out of range\n";
		cerr<<"At this point: "<<__FILE__<<" "<<__LINE__<<endl;
		exit(1);
	}
	if (distance==1) return Nidx[n+index*z(index)];
	if (!(latticeType=="square") || distance>2) {
		cerr<<"Geometry::z: distance>2 or distance==2 && type!=square\n";
		cerr<<"             (not implemented)\n";
		cerr<<"At this point: "<<__FILE__<<" "<<__LINE__<<endl;
		exit(1);
	}
	return Nidx[n+index*z(index)];
}

int Geometry::borderId(int index,int n) const 
{
	if (border.size()<=0) {
		if (!isWarned) {
			cerr<<"WARNING: Geometry::borderId not implemented\n";
			isWarned=true;
		}
		return -1;
	}
	if (index<0 || index>=volume()) {
		cerr<<"Geometry::z: index="<<index<<" out of range\n";
		cerr<<"At this point: "<<__FILE__<<" "<<__LINE__<<endl;
		exit(1);
	}
	return border[n+index*z(index)];
}

std::string Geometry::latticeName() const
{
	if (!isInit) return "NULL"; 
	return latticeType; 
}
		
#ifdef NO_DEBUG	
int main()
{
	int l=2;
	vector<siteIndex> lvector;
	lvector.push_back(l); lvector.push_back(l); lvector.push_back(l); // a 4x4x4 lattice
	Geometry g;
	//g.init("rectangle",lvector);
	//g.init("cubic","2,2,2",1);
	//cerr<<"33+40="<<g.add(33,40)<<endl;
	g.init("fcc","2,2,2",1);
	unsigned int i,j;
	/* for (i=0;i<g.volume();i++) {
		for (j=0;j<g.volume();j++) {
			cerr<<i<<"+"<<j<<"="<<g.add(i,j)<<endl;
		}
	} */
//	cerr<<g.add(8,16)<<endl;
	
	
	return 0;
}
#endif

void vectorSet(vector<double> &vdest,double v)
{
	unsigned int i;
	for (i=0;i<vdest.size();i++) {
		vdest[i] = v;
	}
}
