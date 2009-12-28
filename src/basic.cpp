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


#include "basic.h"
#include "KISS.h"

psimag::KISS rngKiss;


using namespace std;

int mySign(double x)
{
	if (x>0) return 1;
	return -1;
}

void invert(char *temp,int n)
{
	int i;
	char auxtmp[256];
	for (i=0;i<n;i++) {
		auxtmp[i]=temp[n-i-1];
	}		
	
	for (i=0;i<n;i++) {
		temp[i]=auxtmp[i];
	}
}

/*!\brief Transforms an integer into a char *.
	* \param temp : The char * (output)
	* \param n : the number of digits (output)
	* \param x : the integer (input)
	*/
void num2digits(char *temp,int &n,int x)
{
	int r,y,tmp;
	int i=0;
	if (x==0) {
		n=0;
		temp[0]=x+48;
		temp[1]='\0';
		return;
	}
	n = int(log(double(x))/log(10.0));
	y =x;
	while(1) {
		tmp = int(y/10);
		r = y % 10;
		y =tmp;
		temp[i++]=r+48;
		if (y==0) break;
	}
	invert(temp,i);
	temp[i]='\0';
}

/*!\brief Transforms an integer into a string and returns it using num2digits() defined before.
	* \param x : the integer (input)
	*/
string dtos(int x)
{
	char temp[256];
	int n;
	num2digits(temp,n,x);
	string s = temp;
	return s;
}

double fermiZero(double x)
{
	if (x>0) return 0;
	return 1;
}

double fermi(double x)
{
	double res;
	if (x>50) return 0;
	if (x<-50) return 1;
	res = 1.0/(1.0+exp(x));
	return res;
}

double logfermi(double x)
{
	double res;
	if (x>20) return -x;
	if (x<-20) return 0;
	res = -log(1.0+exp(x));
	return res;
}

// Derivative (prime) of Fermi's function
double fermiPrime(double x)
{
	double res;
	res= -fermi(x)*fermi(-x);
	return res;
}		

void tmpValues(double &a,double &b,double &mu,double &beta,int option)
{
	static double a1,b1,mu1,beta1;
	if (option==0) { // set static members
		a1 = a;
		b1 = b;
		mu1 = mu;
		beta1 = beta;
	} else { //retrieve static members
		a = a1;
		b = b1;
		mu = mu1;
		beta = beta1;
	}
}
			 
double maxElement(vector<double> const &v)
{
	vector<double>::const_iterator it = max_element(v.begin(), v.end());
	return *it;
}

double maxElement(vector<tpem_t> const &v)
{
	vector<double> vv(v.size());
	unsigned int i;
	for (i=0;i<v.size();i++) vv[i]=sqrt(real(v[i])*real(v[i])+imag(v[i])*imag(v[i]));
	vector<double>::const_iterator it = max_element(vv.begin(), vv.end());
	return *it;
}

double minElement(vector<double> const &v)
{       
        vector<double>::const_iterator it = min_element(v.begin(), v.end());
        return *it;
}

bool isInVector(int i,vector<int> const &v)
{
        unsigned int j;
        for (j=0;j<v.size();j++) {
                if (v[j]==i) return true;
        }
        return false;
}

void printProgress(int i,int total,int nMarks,char mark,int option)
{
	int every=total/nMarks;
	if (every<=0 || i<=0) return;
	if (i%every ==0) {
		std::cout<<mark;
		std::cout.flush();
	}
}


int parity(int i, int d, int L)
{
	int x=0,y=0,z=0;
	
	switch (d) {
		case 1:
			x=i; y=0; z=0;
			break;
		case 2:
			x=i%L;
			y=(i-x)/L;
			z=0;
			break;
		case 3:
			x=i%L;
			z=(i-x)/L;
			y=z % L;
			z=(z-y)/L;
			break;
	}
	x += (y+z);
	if (x%2==0) return 1;
	else return 0;
}
/*

void myRandomSeed(unsigned  int seed)
{
	std::cout<<"Intializing random number generator with "<<seed;
	std::cout<<" RAND_MAX="<<RAND_MAX<<endl;
	srandom(seed);
}

double myRandom()
{
	//double r_norm = 1.0/(double)0x7fffffff;
	double r_norm = 1.0/RAND_MAX;
	return r_norm*(double)random();
} 
*/

void myRandomSeed(unsigned  int seed)
{
	if (rngKiss.getLength()!=4294967295) {
		cerr<<"SERIOUS PROBLEM: RANDOM NUMBER GENERATOR NOT RELIABLE ON THIS PLATFORM!!\n";
		cerr<<"At this point: "<<__FILE__<<" "<<__LINE__<<endl;
		cerr<<"I will exit now, you should call your customer support 1-800 number.\n";
		exit(1);
	}
	rngKiss.seed(seed);
	
}

double myRandom()
{
	return rngKiss();
}


double vectorSum(vector<double> &v)
{
	//* Should be: foreach(v.begin(),v.end(),Sum); but Sum not defined.
	double s=0.0;
	for (unsigned int i=0;i<v.size();i++) s+= v[i];
	return s;
}

void randomize(vector<double> &v,double v0,double dv)
{
	for (unsigned i=0;i<v.size();i++) {
		v[i] = v0 + dv*(0.5-myRandom());
	}
}

void randomizeBox(vector<double> &v,double v0, double dv)
{
	unsigned int i;
	for (i=0;i<v.size();i++) {
        	if (myRandom()>0.5) {
                	v[i] = v0 + dv; 
        	} else {
                	v[i] = v0 - dv;
        	}
	}
}


void vectorPrint(vector<int> const &v,const std::string& name,ostream &s)
{
	int i,n=v.size();
	//for (i=0;i<n;i++) s<<name<<"["<<i<<"]="<<v[i]<<endl;
	s<<"#"<<name<<endl;
	for (i=0;i<n;i++) s<<i<<" "<<v[i]<<endl;
}

void vectorPrint(vector<double> const  &v,const std::string& name,ostream &s)
{
	int i,n=v.size();
	//for (i=0;i<n;i++) s<<name<<"["<<i<<"]="<<v[i]<<endl;
	s<<"#"<<name<<endl;
	for (i=0;i<n;i++) s<<i<<" "<<v[i]<<endl;
}

void vectorPrint(vector<tpem_t> const &v,const std::string& name,ostream &s)
{
	int i,n=v.size();
	//for (i=0;i<n;i++) s<<name<<"["<<i<<"]="<<v[i]<<endl;
	s<<"#"<<name<<endl;
	for (i=0;i<n;i++) s<<i<<" "<<real(v[i])<<" "<<imag(v[i])<<endl;
}

 
void vectorDivide(vector<double> &v,double value)
{
	if (value==0) return;
	int i,n=v.size();
	for (i=0;i<n;i++) {
		v[i] /= value;
	}
}

void vectorDivide(vector<tpem_t> &v,double value)
{
	if (value==0) return;
	int i,n=v.size();
	for (i=0;i<n;i++) {
		v[i] = v[i]/value;
	}
}

double vectorNorm2(vector<double> const &v)
{
	double ret = 1.0;
	unsigned int i;
	
	for (i=0;i<v.size();i++) {
		if (v[i]==0.0) return 0.0;
		ret *= v[i];
	}
	return ret;
}

template <class T>
void mysplit(string const &s,vector<T> &c,char schar)
{
	int i,l;
	string buffer;
	
	c.clear();
	l=s.length();
	for (i=0;i<l;i++) {
		if (s.at(i)==schar) {
			c.push_back(atof(buffer.c_str()));
			buffer="";
		} else {
			buffer=buffer + s.at(i);
		}
	}
	c.push_back(atof(buffer.c_str()));
	//cerr<<"Split function: tried to split "<<s<<endl;
	//cerr<<"Split function: got "<<c[0]<<" "<<c[1]<<endl;
}

void mychop(string &s)
{
	std::string buffer;
	int i,l=s.length();
	
	for (i=0;i<l-1;i++) {
		buffer = buffer + s.at(i);
	}
	buffer = s;
}

/*!\brief Loads a vector from a file.
	* \param v : The vector that must have been previously allocated (output).
	* \param filename : The name of the file from where to read the vector. It should contain
		a label and then a two column index value per line with the values of
		the vector entries (input).
	* \param label : The label that indicates from where to start to read the file (input).
	* \param level : If there are many labels read the the level-th one.
	*/
int loadVector(vector<double> &v,string const &filename,string &label,int level)
{
	std::ifstream fin(filename.c_str());
	int i,n;
	string tmpLine;
	int maxline=10240;
	char *buffer;
	
	buffer  = new char[maxline];
		
	if (!fin || fin.bad()) {
		cerr<<"Cannot open file "<<filename<<endl;
		cerr<<"AT THIS POINT "<<__FILE__<<" "<<__LINE__<<endl;
		exit(1); // throw an exception instead, caller must catch it
	}
	
	//cerr<<"load vector: I just successfully opened: "<<filename<<endl;
	for (i=0;i<level;i++) { // .sav will have 2 data sets, read the level-th data set
		while(!fin.eof()) {
			fin.getline(buffer,maxline);
			tmpLine=string(buffer);
			if (fin.fail()) {
				cerr<<"loadVector: Premature end of file="<<filename<<" no label="<<label<<" found.\n";
				cerr<<"AT THIS POINT "<<__FILE__<<" "<<__LINE__<<endl;
				exit(1);
			}
			//cerr<<"Load vector, before chopping "<<tmpLine<<"*"<<endl;
			mychop(tmpLine);
			//cerr<<"Load vector, after chopping "<<tmpLine<<"*"<<endl;
			if (tmpLine==label) {
				//cerr<<"I read "<<tmpLine<<" which IS equal to "<<label<<endl;
				break;
			}
			//cerr<<"I read "<<tmpLine<<" which is not equal to "<<label<<endl;
		}
	}
	
	if (!(tmpLine==label)) {
		cerr<<"loadVector: Problem trying to read file "<<filename<<endl;
		cerr<<"Label="<<label<<" not found\n";
		cerr<<"AT THIS POINT "<<__FILE__<<" "<<__LINE__<<endl;
		exit(1);
	}
	n = v.size();
	vector<double> c;
	
	for (i=0;i<n;i++) {
		fin.getline(buffer,maxline);
		tmpLine=string(buffer);
		mysplit(tmpLine,c,' ');
		if (i!=c[0]) {
			cerr<<"loadVector: Mismatch while reading label="<<label<<" in file "<<filename;
			cerr<<" "<<i<<" not equal "<<c[0]<<endl;
			cerr<<"AT THIS POINT "<<__FILE__<<" "<<__LINE__<<endl;
			exit(1);// throw an exception instead, caller must catch it
		} 
		v[i]=c[1];
	}
	fin.close();
	delete [] buffer;
	return 0;
}
	
void randomModulus(vector<int> &v,int conc,int n)
{
	unsigned int i,j;
	v.assign(v.size(),0);
	for (i=0;i<conc;i++) {
		j=0;
		while (v[j]==1) j=int(n*myRandom());		
		v[j]=1;
	}
}
template void mysplit(std::string const &s,vector<int> &v,char schar);
template void mysplit(std::string const &s,vector<double> &v,char schar);

