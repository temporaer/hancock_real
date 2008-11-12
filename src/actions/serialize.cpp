#include <dlfcn.h>
#include <boost/shared_ptr.hpp>
#include <factory/factory.h>
#include <cstdio>

#include <configuration.hpp>
#include <sdp_wrapper.hpp>
#include <seriation_gen.hpp>
#include <adjmat_gen.hpp>

#include "sdp_simple.hpp"
#include <nana.h>

using namespace std;
using namespace boost;

void SDPSimple::operator()()
{
	string adjmat_gen_name               = "RandomAdjMatGen";
	auto_ptr<AdjMatGen> adjmat_gen       = genericFactory<AdjMatGen>::instance().create(adjmat_gen_name);
	if(!adjmat_gen.get())
		throw logic_error(string("Supplied AdjMatGen `") + adjmat_gen_name + "' does not exist");
	adjmat_gen->configure();

	shared_ptr<AdjMatGen::AdjMatT> adjmat_ptr((*adjmat_gen)());

	string seriation_gen_name            = "SDPSeriationGen";
	auto_ptr<SeriationGen> seriation_gen = genericFactory<SeriationGen>::instance().create(seriation_gen_name);

	seriation_gen->configure();

	SeriationGen::SeriationT randwalk;
	randwalk = (*seriation_gen)(adjmat_ptr);
}


namespace{
	registerInFactory<Action, SDPSimple> registerBase("SDPSimple");
}
