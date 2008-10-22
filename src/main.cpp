/*       Created   :  10/03/2008 08:22:01 PM
 *       Last Change: Wed Oct 22 09:00 PM 2008 CEST
 */

#include <dlfcn.h>
#include <string>
#include <iostream>
#include <factory/factory.h>
#include "configuration.hpp"

#include <sdp_seriation_gen.hpp>
#include <sdp_wrapper.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <nana.h>
using namespace boost::numeric::ublas;
using namespace std;

int main(int argc, char* argv[]){
	dlopen("sdp_wrappers/libsdp_wrappers.so",RTLD_LAZY);
	
	gCfg().parsecfg(argc,argv);

	typedef matrix<double,column_major> AdjMatCT;
	AdjMatCT adjmat(6,6);
	for(int i=0;i<6;i++)
        for(int j=0;j<6;j++)
            adjmat(i,j) = ((i+j)%2>0)?1:0;

	SDPSeriationGen walkgen;
	string sdp_wrapper_name          = gCfg().getString("ser-sdp-wrapper");
	auto_ptr<SDPWrapper> sdp_wrapper = genericFactory<SDPWrapper>::instance().create(sdp_wrapper_name);
	walkgen.setSDPWrapper( sdp_wrapper );

	SeriationGen::SeriationT randwalk  = walkgen(adjmat);
}
