/*       Created   :  10/07/2008 12:08:01 PM
 *       Last Change: Tue Oct 07 09:00 PM 2008 CEST
 */

#include <configuration.hpp>
#include <boost/program_options.hpp>
#include <factory/factory.h>

struct SDPWrapperCfg{
	SDPWrapperCfg();
};

using namespace boost::program_options;
using namespace std;
SDPWrapperCfg::SDPWrapperCfg(){
	options_description OD;
#ifdef HAVE_DSDP
	options_description dsdp("DSDP Options");
	dsdp.add_options()
		("f",value<float>(),"f, a float.")
		;
	OD.add(dsdp);
#endif
#ifdef HAVE_CSDP
	options_description csdp("CSDP Options");
	csdp.add_options()
		("g",value<float>(),"g, a float.")
		;
	OD.add(csdp);
#endif
	gCfg().addModuleOptions(OD);
}
namespace{
	SDPWrapperCfg _sdpwrappercfg;
}
