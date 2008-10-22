/*       Created   :  10/06/2008 12:52:01 AM
 *       Last Change: Wed Oct 22 08:00 PM 2008 CEST
 */

#include <fstream>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <nana.h>
#include <configuration.hpp>
#include <sdp_prob.hpp>
#include <dsdp_wrapper.hpp>
#include <dsdp/dsdp5.h>
#include <factory/factory.h>
using namespace std;
using namespace boost::numeric;

#ifndef DSDP_BINARY
#  error "You need to define DSDP_BINARY, pointing me to where dsdp resides"
#endif

struct DSDPWrapper::Impl{
	typedef DSDPWrapper::AnswerT AnswerT;
	AnswerT operator()(const SDPProb&);
	void writeSDPASparseInputFile(const SDPProb& p, const char* fn);
};

/*!
* helper function to print a matrix in SDPA Sparse format
*/
template<class T>
void sdpaSparsePrintMat(ofstream& o,int n, const T& m)
{
	for(unsigned int i=0;i<m.size1();i++)
		for(unsigned int j=i;j<m.size2();j++)
			o << n << " 1 " << (i+1) << " " << (j+1) << " " << m(i,j)<< endl;
}


void DSDPWrapper::Impl::writeSDPASparseInputFile(const SDPProb& p, const char* fn)
{
	ofstream o(fn);
	o << "\"SDPASparse Sample Problem\"" << endl
	<< p.F.size()                        << endl  // num of constraint matr
	<< "1"                               << endl  // number of blocks in SDP
	<< p.C.size1()                       << endl; // size of the blocks from prev line (one line per block)
	if (p.b.size() > 0)                            // b, the righthandside vec, all on one line
		o << p.b(0);
	for (unsigned int i=1;i<p.b.size();i++)
		o << " " << p.b(i);
	o << endl;

	sdpaSparsePrintMat(o,0, p.C);
	// TODO: Identity-Matrix could be handled separatly (just specify diag)
	for (unsigned int f=0;f<p.F.size();f++) {
	    const SDPProb::MatT& m = p.F[f];
	    sdpaSparsePrintMat(o,f+1,m);
	}
}

bool readDSDPOutputFile(const char* out, DSDPWrapper::AnswerT& ret){
	ifstream is(out);
	string line;
	getline(is,line);
	boost::trim(line);

	vector<string> strvec;
	boost::split( strvec, line, boost::is_any_of(" ") );

	// convert strings to doubles
	int i=0;
	ret = DSDPWrapper::AnswerT(strvec.size());
	for(vector<string>::iterator it = strvec.begin();it!=strvec.end();it++,i++){
		L(" --> %s \n", it->c_str());
		ret(i) = boost::lexical_cast<double>(*it);
	}
	return true;
}

void runDSDP(const char* in, const char* out)
{
	char cmd[255];
	sprintf(cmd,"%s %s -save %s 2>&1 > /dev/null",DSDP_BINARY,in,out);
	L("Executing cmd: %s",cmd);
	system(cmd);
}


DSDPWrapper::AnswerT DSDPWrapper::Impl::operator()(const SDPProb& p){
	L("DSDPWrapper::operator()\n");
	const char* fn_in  = "/tmp/x.dat";
	const char* fn_out = "/tmp/x.out";
	writeSDPASparseInputFile(p,fn_in);
	runDSDP(fn_in,fn_out);
	DSDPWrapper::AnswerT ret;
	readDSDPOutputFile(fn_out,ret);
	return ret;
}


// Wrappers
DSDPWrapper::DSDPWrapper()
	: mImpl(new Impl)
{
}
DSDPWrapper::~DSDPWrapper()
{
    L("Destroying DSDPWrapper");
}
DSDPWrapper::AnswerT DSDPWrapper::operator()(const SDPProb& prob)
{
	return (*mImpl)(prob);
}

namespace{
	registerInFactory<SDPWrapper, DSDPWrapper>  registerBase("DSDPWrapper");
}
