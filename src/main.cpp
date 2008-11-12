/*       Created   :  10/03/2008 08:22:01 PM
 *       Last Change: Wed Nov 12 07:00 PM 2008 CET
 */

#include <dlfcn.h>
#include <string>
#include <iostream>
#include <fstream>
#include <factory/factory.h>
#include "configuration.hpp"
#include "action.hpp"

#include <sdp_seriation_gen.hpp>
#include <sdp_wrapper.hpp>

#include <boost/shared_ptr.hpp>
#include <nana.h>
#undef C
using namespace boost::numeric::ublas;
using namespace boost;
using namespace std;


int main(int argc, char* argv[]){
	dlopen("actions/libactions.so",RTLD_LAZY);

	gCfg().parsecfg(argc,argv);

	string action_name = gCfg().getString("action");
	auto_ptr<Action> action = genericFactory<Action>::instance().create(action_name);
	if(!action.get())
		throw logic_error(string("Supplied action `") + action_name + "' does not exist");
	(*action)();
}
