/*       Created   :  10/06/2008 12:52:01 AM
 *       Last Change: Tue Oct 07 08:00 PM 2008 CEST
 */

#include <nana.h>
#include <csdp_wrapper.hpp>
#include <factory/factory.h>
CSDPWrapper::AnswerT CSDPWrapper::operator()(const SDPProb&)
{
	L("CSDPWrapper::operator()\n");
	AnswerT ret;
	return ret;
}
 CSDPWrapper::~CSDPWrapper(){
L("Destroying CSDPWrapper\n");
 }

namespace{
	registerInFactory<SDPWrapper, CSDPWrapper>  registerBase("CSDPWrapper");
}
