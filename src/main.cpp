/*       Created   :  10/03/2008 08:22:01 PM
 *       Last Change: Thu Oct 23 10:00 AM 2008 CEST
 */

#include <dlfcn.h>
#include <string>
#include <iostream>
#include <factory/factory.h>
#include "configuration.hpp"

#include <sdp_seriation_gen.hpp>
#include <sdp_wrapper.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/shared_ptr.hpp>
#include <nana.h>
using namespace boost::numeric::ublas;
using namespace boost;
using namespace std;

int main(int argc, char* argv[]){
	dlopen("sdp_wrappers/libsdp_wrappers.so",RTLD_LAZY);
	
	gCfg().parsecfg(argc,argv);

	int n=100;
	typedef matrix<double,column_major> AdjMatCT;
	shared_ptr<AdjMatCT> adjmat_ptr(new AdjMatCT(n,n));
	AdjMatCT& adjmat = *adjmat_ptr;
	for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
            adjmat(i,j) = ((i+j)%2>0)?1:0.0;

	SDPSeriationGen walkgen;
	string sdp_wrapper_name          = gCfg().getString("ser-sdp-wrapper");
	auto_ptr<SDPWrapper> sdp_wrapper = genericFactory<SDPWrapper>::instance().create(sdp_wrapper_name);
	sdp_wrapper->configure();
	walkgen.setSDPWrapper( sdp_wrapper );

	SeriationGen::SeriationT randwalk  = walkgen(adjmat_ptr);
}
