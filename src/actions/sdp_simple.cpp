#include <dlfcn.h>
#include <boost/shared_ptr.hpp>
#include <factory/factory.h>
#include <cstdio>

#include <configuration.hpp>
#include <sdp_wrapper.hpp>
#include <seriation_gen.hpp>
#include <random_adjmat_gen.hpp>

#include "sdp_simple.hpp"
#include <nana.h>

using namespace std;
using namespace boost;

void SDPSimple::operator()()
{
	dlopen("sdp_wrappers/libsdp_wrappers.so",RTLD_LAZY);
	RandomAdjMatGen matgen;
	matgen.configure();
	shared_ptr<RandomAdjMatGen::AdjMatT> adjmat_ptr(matgen());

	string seriation_gen_name            = "SDPSeriationGen";
	auto_ptr<SeriationGen> seriation_gen = genericFactory<SeriationGen>::instance().create(seriation_gen_name);

	seriation_gen->configure();

	SeriationGen::SeriationT randwalk;
	randwalk = (*seriation_gen)(adjmat_ptr);
}


namespace{
	registerInFactory<Action, SDPSimple> registerBase("SDPSimple");
}
