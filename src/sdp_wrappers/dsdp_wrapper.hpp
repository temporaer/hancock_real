/*       Created   :  10/06/2008 12:57:46 AM
 *       Last Change: Mon Oct 06 12:00 AM 2008 CEST
 */

#ifndef __DSDP_WRAPPER_HPP__
#define __DSDP_WRAPPER_HPP__
#include <sdp_wrapper.hpp>
class DSDPWrapper: public SDPWrapper
{
	virtual AnswerT operator()(const SDPProb&);
};
#endif

