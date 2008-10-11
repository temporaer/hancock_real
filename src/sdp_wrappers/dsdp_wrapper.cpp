/*       Created   :  10/06/2008 12:52:01 AM
 *       Last Change: Tue Oct 07 08:00 PM 2008 CEST
 */

#include <nana.h>
#include <configuration.hpp>
#include <dsdp_wrapper.hpp>
#include <dsdp/dsdp5.h>
#include <factory/factory.h>

struct DSDPWrapper::Impl{
	typedef DSDPWrapper::AnswerT AnswerT;
	AnswerT operator()(const SDPProb&);
};


DSDPWrapper::AnswerT DSDPWrapper::Impl::operator()(const SDPProb&prob){
	L("DSDPWrapper::operator()\n");
	// TODO: wrap dsdp.
	AnswerT ret;
	return ret;
}


// Wrappers
DSDPWrapper::DSDPWrapper()
	: mImpl(new Impl)
{
}
DSDPWrapper::~DSDPWrapper()
{
    L("Destroying DSDPWrapper");
}
DSDPWrapper::AnswerT DSDPWrapper::operator()(const SDPProb& prob)
{
	return (*mImpl)(prob);
}

namespace{
	registerInFactory<SDPWrapper, DSDPWrapper>  registerBase("DSDPWrapper");
}
