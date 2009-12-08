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




#ifndef GEOMETRY_H_
#define GEOMETRY_H_
#include <vector>
#include <iostream>
#include <string>
#include <cmath>
#include "basic.h"

using std::vector;
using std::string;
using std::endl;

extern void vectorSet(vector<double> &vdest,double v);

inline void vectorPrint(vector<double> const &v,string const &label)
{
	unsigned int i;
	for (i=0;i<v.size();i++) std::cerr<<label<<"["<<i<<"]="<<v[i]<<endl;
}

inline vector<double> operator*(vector<double> const &a,double v) 
{
	vector<double> c(a.size());
	unsigned int i;
	for (i=0;i<a.size();i++) c[i] = a[i]*v;
	return c;
}

inline vector<double> operator*(double v,vector<double> const &a) 
{
	vector<double> c(a.size());
	unsigned int i;
	for (i=0;i<a.size();i++) c[i] = a[i]*v;
	return c;
}

inline vector<double> operator+(vector<double> const &a,vector<double> const &b) 
{
	vector<double> c(a.size());
	unsigned int i;
	for (i=0;i<a.size();i++) c[i] = a[i]+b[i];
	return c;
}
inline vector<double> operator-(vector<double> const &a,vector<double> const &b) 
{
	vector<double> c(a.size());
	unsigned int i;
	for (i=0;i<a.size();i++) c[i] = a[i]-b[i];
	return c;
}
typedef  int siteIndex;

/*! \brief Geometry class to handle neighbors
 * Replaces old Lattice class
 */
class Geometry {

public:
	/** Default constructor so it can be called early  */
	Geometry() { isInit=false; errorPrinted=0; isWarned=false; } 
	~Geometry() { isInit=false; } 
	void init(string const &s1,string const &s2,int verbose=0);
	/** \brief Returns the number of neighbors of site index at distance */
	int z(int index=0,int distance=1) const;
	/** Returns the n-th neighbor of site index at distance */
	int neighbor(int index,int n,int distance=1) const ;
	/** Length of the lattice along dimension d */
	int length(int d) const { return Length[d]; }
	/** Number of lattice sites  */
	int volume() const { return Nn.size(); }
	/** Prints the neighbors of all sites. Usefull for debugging */
	void nicePrint(std::ostream &s,int option=0) const;
	/** \brief Returns the dimension of the lattice */
	int dim() const;
	/** \brief Returns the length of the lattice. Implemented only for lattices with all sides equal. */
	int length() const;
	/** \brief Lattice Addition.
	 
	 * Converts ind and ind2 to vectors on the lattice, adds them, and returns the corresponding
	 * index.  Implemented only for lattices with all sides equal. */
	int add(siteIndex ind,siteIndex ind2) const;
	/** \brief Scalar Product
	
	 * Converts i and j to vectors on the lattice, performs the scalar product and returns the 
	 * corresponding index. Implemented only for lattices with all sides equal. */
	double scalarProd(siteIndex i,siteIndex j) const;
	
	/** \brief Border function
	
	* If the j-th neighbor of i does not cross a lattice boundary it returns -1.
	* Otherwise it returns the direction across which a lattice boundary is crossed */
	int borderId(siteIndex i,siteIndex j) const;
	
	/** \brief Returns the type of Lattice as initialized with init*/
	std::string latticeName()  const;
	
	/** \brief Returns true if the lattice is "1d", "square" or "cubic" and false otherwise */
	bool isCubicType(string const &s) const;
	bool isCubicType() const;
	void direction(vector<double> &r,siteIndex site1,siteIndex site2) const;
	int scalarDirection(siteIndex site1,siteIndex site2,int number= -1) const;
	
private:
	struct g_aux_Zb_data
	{
		double minDistance;
		int length;
		int volume,nbasis;
	};

	mutable int errorPrinted;
	string latticeType; 
	vector<siteIndex> Length; /**< stores the lattice sides */
	vector<siteIndex> Nidx; /**<stores the neighbors at distance =1 */
	vector<siteIndex> Nn; /**<stores the number of neighbors at distance =1*/
	vector<siteIndex> Nidx2;
	/**<stores the neighbors at distance =1 (not available for all lattices)  */
	vector<siteIndex> Nn2; 
	/**<stores the number of neighbors at distance =2 (not available for all lattices) */
	vector<int> border;
	/**<which border is crossed if any (not available for all lattices) */
	bool isInit; /**< indicates if object has been initialized. */
	mutable bool isWarned;
	
	/** Initializer for prism/rectangle lattices  */
	void init(string const &type,vector<siteIndex> const &Lvector);
	/** Initializer for cubic/square/1d lattices  */
	void init(string const &type,int l);
	
	bool g_pbc(int &x, int L) const;
	int g_index(int &x,int &y,int &z,int L);
	int g_index(int &x,int &y,int &z,vector<siteIndex> const &lvector);
	void g_geometryBcc(vector<siteIndex> &nn,vector<siteIndex> &nidx,int l);
	void g_add(int &xx,int &yy, int &zz,vector<int> const &v) const;
	void g_add(vector<double> &dest,vector<double> const &src) const;
	void g_auxFcc(vector<siteIndex> &nidx,vector< vector<int> > const &v,
		int x,int y,int z,int l,int bshift);
	void g_auxFcc(vector< vector<int> > &v,vector< vector<int> > const &w);
	void g_geometryFcc(vector<siteIndex> &nn,vector<siteIndex> &nidx,int l);
	void g_geometryTriangular(vector<siteIndex> &nn,vector<siteIndex> &nidx,int l);
	void g_geometryHoneycomb(vector<siteIndex> &nn,vector<siteIndex> &nidx,vector<int> const &lvector);
	void g_geometryPrism(vector<siteIndex> &nn,vector<siteIndex> &nidx,vector<siteIndex> const &lvector);
	void g_geometryCubic(vector<siteIndex> &nn,vector<siteIndex> &nidx,int l);
	void g_geometryRectangle(vector<siteIndex> &nn,vector<siteIndex> &nidx,
		vector<siteIndex> const &lvector);
	void g_geometryZb(vector<siteIndex> &nn,vector<siteIndex> &nidx,int l);
	void g_geometryOneD(vector<siteIndex> &nn,vector<siteIndex> &nidx,int l);
	void g_geometrySquare(vector<siteIndex> &nn,vector<siteIndex> &nidx,int l);
	void geometry(vector<siteIndex> &nn,vector<siteIndex> &nidx,int l,string const &type);
	void geometry(vector<siteIndex> &nn,vector<siteIndex> &nidx,
		vector<siteIndex> const &lvector,string const &type);
	void g_geometryDistance2(vector<siteIndex> &nn2,vector<siteIndex> &nidx2,
		vector<siteIndex> const &lvector);
	/** \brief Functions valid only for lattices with all lattices equal */
	siteIndex coor2Index(vector<siteIndex> &v,string const &s) const;
	//template<class T>
	//void index2Coor(vector<T> &v,siteIndex i,std::string const &lt) const;
	void index2Coor(vector<double> &v,siteIndex i,std::string const &lt) const;
	void index2Coor(vector<siteIndex> &v,siteIndex i,std::string const &lt) const;
	
	bool allSidesEqual() const;
	
	double g_norm(vector<double> const &r);
	double g_auxZb_dist(siteIndex site1, siteIndex site2,g_aux_Zb_data const &data);
	void g_auxZb_initBasis(vector<vector<double> > &basisVector) const;
	double g_auxZb_distance(vector<double> const &r1,vector<double> const &r2,g_aux_Zb_data const &data);
	void g_auxZb_Index2Coor(vector<double> &r,int site,int volume) const;
	void g_auxZb_findNn(vector<siteIndex> &v,int site,g_aux_Zb_data const &data);
	
	double g_auxFcc_dist(siteIndex site1, siteIndex site2,g_aux_Zb_data const &data);
	void g_auxFcc_initBasis(vector<vector<double> > &basisVector) const;
	double g_auxFcc_distance(vector<double> const &r1,vector<double> const &r2,g_aux_Zb_data const &data);
	void g_auxFcc_Index2Coor(vector<double> &r,int site,int vol) const;
	void g_auxFcc_findNn(vector<siteIndex> &v,int site,g_aux_Zb_data const &data);
	void g_geometryFccDirection(vector<double> &r,siteIndex site1, siteIndex site2) const;
	void g_geometryZbDirection(vector<double> &r,siteIndex site1, siteIndex site2) const;
	
	int add_fcc_aux(vector<double> const &b1,vector<double> const &b2) const;
	void add_fcc(vector<double> const &r1,vector<double> const &b1,vector<double> const &r2,vector<double> const &b2,
vector<double> &rsum,vector<double> &bsum) const;
	int add_fcc(siteIndex ind,siteIndex ind2) const;
};
	

#endif
