/*       Created   :  10/06/2008 12:52:01 AM
 *       Last Change: Sun Oct 12 10:00 PM 2008 CEST
 */
#include <fstream>
#include <cstdlib>
#include <nana.h>
#include <sdp_prob.hpp>
#include <sdpa_wrapper.hpp>
#include <factory/factory.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <configuration.hpp>

#ifndef SDPA_BINARY
#  error "You need to define SDPA_BINARY, pointing me to where sdpa resides"
#endif

using namespace std;

template<class T>
void sdpaPrintMat(ofstream& o,const T& m)
{
	SDPProb::MatT::const_iterator1 row=m.begin1(),row_begin=m.begin1(),row_end=m.end1();
	SDPProb::MatT::const_iterator2 col;
	o << "{";
	for (row = m.begin1(); row != row_end; row++) {
		if (row!=row_begin) o << ", ";
		o << "{";
		copy(row.begin(),row.end(),ostream_iterator<double>(o,", "));
		o << "}";
	}
	o << "}"<<endl;
}

void SDPAWrapper::writeSDPAInputFile(const SDPProb& p, const char* fn)
{
	ofstream o(fn);
	o   << "\"SDPA-Wrapper Problem\""         << endl
	<< " "<<p.C.size1()<<" = mDIM\n"      << endl
	<< " 1 = nBLOCK"                      << endl
	<< " "<<p.C.size1()<<" = bLOCKsTRUCT" << endl
	<< "{";
	if (p.b.size() > 0)
		o << p.b(0);
	for (unsigned int i=1;i<p.b.size();i++)
		o << ", " << p.b(i);
	o   << "}"                                << endl;
	sdpaPrintMat(o,p.C);
	// TODO: Identity-Matrix could be handled separatly (just specify diag)
	for (unsigned int f=0;f<p.F.size();f++) {
	    const SDPProb::MatT& m = p.F[f];
	    sdpaPrintMat(o,m);
	}
}
bool SDPAWrapper::readSDPAOutputFile(const char*out, AnswerT& ret)
{
	ifstream is(out);
	string line;
	bool found = false;
	while( !found && getline(is,line) ){
		string start ( line.substr(0,6) );
		if( start == "xVec =" )
			found = true;
	}
	if(found == false){
		IP(found,"Could not parse SDPA output, xVec not found");
		return false;
	}
	getline(is,line); // get xVec
	line = line.substr(1,line.size()-2);


	// extract numbers as strings
	vector<string> strvec;
	boost::split( strvec, line, boost::is_any_of(",") );

	// convert strings to doubles
	int i=0;
	ret = AnswerT(strvec.size());
	for(vector<string>::iterator it = strvec.begin();it!=strvec.end();it++,i++){
		ret(i) = boost::lexical_cast<double>(*it);
	}

	return true;
}

void SDPAWrapper::runSDPA(const char* in, const char* out, const char* param)
{
	char cmd[255];
	sprintf(cmd,"%s -p %s -dd %s -o %s 2>&1 > /dev/null",SDPA_BINARY,param,in,out);
	system(cmd);
}

SDPAWrapper::AnswerT SDPAWrapper::operator()(const SDPProb&p)
{
	L("SDPAWrapper::operator()\n");
	const char* fn_in  = "/tmp/x.dat";
	const char* fn_out = "/tmp/x.out";
	const char* fn_par = gCfg().getString("sdpa-param-file").c_str();
	writeSDPAInputFile(p,fn_in);
	runSDPA(fn_in,fn_out,fn_par);
	AnswerT ret;
	readSDPAOutputFile(fn_out,ret);
	return ret;
}
SDPAWrapper::~SDPAWrapper()
{
	L("Destroying SDPAWrapper\n");
}

namespace
{
registerInFactory<SDPWrapper, SDPAWrapper>  registerBase("SDPAWrapper");
}
