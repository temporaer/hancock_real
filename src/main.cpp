/*       Created   :  10/03/2008 08:22:01 PM
 *       Last Change: Tue Oct 28 01:00 AM 2008 CET
 */

#include <dlfcn.h>
#include <string>
#include <iostream>
#include <fstream>
#include <factory/factory.h>
#include "configuration.hpp"

#include <sdp_seriation_gen.hpp>
#include <sdp_wrapper.hpp>

#include <random_adjmat_gen.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/shared_ptr.hpp>
#include <nana.h>
using namespace boost::numeric::ublas;
using namespace boost;
using namespace std;

template<class T,class U>
void writeToFile(
		const char* fn,
		const T&    adj,
		const U&    path
		){
	ofstream o(fn);
	o << "A = [ " ;
	for(int row=0;row<adj.size2();row++){
		for(int col=0;col<adj.size1();col++){
			o << adj(row,col) << " ";
		}
		if(row != adj.size2()-1)
			o << " ; ";
	}
	o << " ];"<<endl;

	o << "path = [ ";
	for(int i=0;i<path.size();i++)
	{
		o << path[i].first << " ";
	}
	o << path.back().second << " ] ; "<< endl;
}

int main(int argc, char* argv[]){
	dlopen("sdp_wrappers/libsdp_wrappers.so",RTLD_LAZY);
	
	gCfg().parsecfg(argc,argv);

	RandomAdjMatGen matgen;
	matgen.configure();
	shared_ptr<RandomAdjMatGen::AdjMatT> adjmat_ptr(matgen());

	SDPSeriationGen walkgen;
	string sdp_wrapper_name          = gCfg().getString("ser-gen.sdp-wrapper");
	auto_ptr<SDPWrapper> sdp_wrapper = genericFactory<SDPWrapper>::instance().create(sdp_wrapper_name);
	sdp_wrapper->configure();
	walkgen.setSDPWrapper( sdp_wrapper );

	SeriationGen::SeriationT randwalk  = walkgen(adjmat_ptr);

	const char* out_file = gCfg().getString("output").c_str();
	writeToFile(out_file,*adjmat_ptr,randwalk);
}
