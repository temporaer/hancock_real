/*       Created   :  10/07/2008 08:30:06 PM
 *       Last Change: Tue Oct 07 08:00 PM 2008 CEST
 */

#include <exception>
#include <nana.h>
#include <sdp_wrapper.hpp>
#include <factory/factory.h>

using namespace std;

SDPWrapper::AnswerT SDPWrapper::operator()(const SDPProb&)
{
	throw logic_error("Called SDPWrapper() w/o subclassing");
}
SDPWrapper::~SDPWrapper(){
    L("Destroying SDPWrapper\n");
}


namespace{
	registerInFactory<SDPWrapper, SDPWrapper> registerBase("SDPWrapper");
}
