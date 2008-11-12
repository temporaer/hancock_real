#include <stdexcept> 
#include <factory/factory.h>
#include "seriation_gen.hpp"
using namespace std;
void SeriationGen::configure()
{
}
SeriationGen::SeriationT SeriationGen::operator()(boost::shared_ptr<AdjMatT>)
{
	throw logic_error("Called SeriationGen() w/o subclassing");
}


namespace{
	registerInFactory<SeriationGen, SeriationGen> registerBase("SeriationGen");
}
