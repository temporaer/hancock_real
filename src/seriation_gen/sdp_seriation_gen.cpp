/*       Created   :  10/05/2008 10:25:05 PM
 *       Last Change: Tue Oct 07 09:00 PM 2008 CEST
 */
#include <nana.h>
#include <sdp_seriation_gen.hpp>
#include <sdp_seriation_prob_gen.hpp>
#include <sdp_wrapper.hpp>
#include <factory/factory.h>

using namespace boost::numeric::ublas;
using namespace std;

// Implementation of SeriationGen

struct SDPSeriationGen::Impl{
    std::auto_ptr<SDPWrapper> mSDPWrapper;
    typedef SeriationGen::AdjMatT AdjMatT;
    SeriationGen::SeriationT operator()(const AdjMatT& adj);
    void setSDPWrapper(std::auto_ptr<SDPWrapper> w);
    ~Impl();
};
SDPSeriationGen::Impl::~Impl(){
}

SeriationGen::SeriationT SDPSeriationGen::Impl::operator()(const AdjMatT& adj)
{
	I(adj.size1() == adj.size2());
	I(mSDPWrapper.get() != NULL);

	// Generate the SDP-Problem
	SDPProb prob;
	SDPSeriationProbGen sdp_probgen(adj);
	sdp_probgen(prob);

	// Solve the SDP-Problem using an SDP Wrapper
	(*mSDPWrapper)(prob);

	SeriationT ret;
	return ret;
}


void SDPSeriationGen::Impl::setSDPWrapper(std::auto_ptr<SDPWrapper> w){
  mSDPWrapper = w;
}




// Wrap SeriationGen::Impl

SDPSeriationGen::SDPSeriationGen()
:mImpl(new Impl)
{

}


SeriationGen::SeriationT SDPSeriationGen::operator()(const AdjMatT& adj)
{
    return (*mImpl)(adj);
}

void SDPSeriationGen::setSDPWrapper(std::auto_ptr<SDPWrapper> w){
    mImpl->setSDPWrapper(w);
}

SDPSeriationGen::~SDPSeriationGen()
{
L("Destroying SDPSeriationGen\n");
}

