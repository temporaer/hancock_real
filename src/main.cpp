/*       Created   :  10/03/2008 08:22:01 PM
 *       Last Change: Mon Oct 27 11:00 PM 2008 CET
 */

#include <dlfcn.h>
#include <string>
#include <iostream>
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
}
