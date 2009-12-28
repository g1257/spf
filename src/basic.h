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




#ifndef _BASIC_H_
#define _BASIC_H_


#ifdef MODEL_KONDO_INF_ONEBAND_PHONONS
#define MODEL_KONDO_PHONONS
#endif

#ifdef MODEL_KONDO_INF_ONEBAND_PHONONS_EX
#define MODEL_KONDO_INF_ONEBAND_PHONONS
#define MODEL_KONDO_PHONONS_EX
#define MODEL_KONDO_PHONONS
#endif

#ifdef MODEL_KONDO_FINITE_PHONONS_EX
#define MODEL_KONDO_FINITE
#define MODEL_KONDO_PHONONS_EX
#define MODEL_KONDO_PHONONS
#endif


#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <string>
#include <sys/types.h>
#include <sys/times.h>
#include <unistd.h>
#include <sys/utsname.h>
#include <cctype>
#include <ctime>
#include <algorithm>
#include <numeric>
#include <vector>
#include <string>
#include "tpemplusTypes.h"
#include "basic_types.h"

using std::vector;


/* File basic.h */
double fermiZero(double x);
double fermi(double x);
double logfermi(double x);
void tmpValues(double &a,double &b,double &mu,double &beta,int option);
// Derivative (prime) of Fermi's function
double fermiPrime(double x);
double maxElement(vector<double> const &v);
double maxElement(vector<tpem_t> const &v);
double minElement(vector<double> const &v);
void printProgress(int i,int total,int nMarks,char mark,int option);
tpem_t complexConvert(MatType const &c);
int parity(int i, int d, int L);
double myRandom();
std::string dtos(int i);
void myRandomSeed(unsigned  int seed);
void randomize(vector<double> &v,double v2,double dv);
void randomizeBox(vector<double> &v,double v0,double dv);
int loadVector(vector<double> &v,std::string const &filename,std::string &label,int
level);
void mychop(std::string &s);
double vectorNorm2(vector<double> const &v);
void vectorDivide(vector<MatType> &v,double value);
void vectorDivide(vector<tpem_t> &v,double value);
void vectorDivide(vector<double> &v,double value);
void vectorPrint(vector<MatType> const &v,const std::string &name,std::ostream &s);
void vectorPrint(vector<tpem_t> const &v,const std::string &name,std::ostream &s);
void vectorPrint(vector<double> const  &v,const std::string &name,std::ostream &s);
void vectorPrint(vector<int> const &v,const std::string& name,std::ostream &s);
double vectorSum(vector<double> &v);
void randomModulus(vector<int> &v,int conc,int n);
bool isInVector(int i,vector<int> const &v);
int mySign(double);
template<class T>
 void mysplit(std::string const &s,vector<T> &v,char schar);

template <class T>
T square(T x)
{
	return (x*x);
}

template<class T>
std::string ttos(T t)
{
	std::stringstream ss;
	std::string str;
	ss<<t;
	ss>>str;
	return str;
}



#endif
