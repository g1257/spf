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





#include "argument.h"

void remove_blanks(std::string &line)
{
	int i,l=line.length();
	std::string buffer ="";
	int blanks=0;
	
	for (i=0;i<l;i++) {
		if (line[i]!=' ' && line[i]!='\t') break;
	}
	while (i<l) {
		if (blanks<2) buffer = buffer + line[i];
		if (line[i]==' ' || line[i]=='\t') {
			blanks ++;
		} else {
			blanks=0;
		}
		i++;
	}
	line = buffer;
}

int procLine(std::string const &line,std::ostream &s)
{
	int i,l=line.length();
	std::string buffer ="";
	bool valueflag=false;
	
	if (line[0]=='#') return 0;
	
	for (i=0;i<l;i++) {
		if (valueflag) {
			buffer = buffer + line[i];
		} 
		if (line[i]==' ' || line[i]=='\t') {
			valueflag=true;
		}
	}
	
	//remove_blanks(buffer);
	s<<buffer<<std::endl;
	return 0;
}
		
int procFile(char const *filename,char const *outfile)
{
	std::ifstream fin(filename);
	int i;
	std::string line;
	std::ofstream fout(outfile);
	
	if (!fin or !fin.good()) {
		std::cerr<<"Cannot open/read file "<<filename<<std::endl;
		return 2;
	}
	while(!fin.eof()) {
		getline(fin,line);
		if ((i=procLine(line,fout))!=0) {
			return i;
		}
	}
	fin.close();
	fout.close();
	return 0;
}
