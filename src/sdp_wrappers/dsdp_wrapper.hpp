/*       Created   :  10/06/2008 12:57:46 AM
 *       Last Change: Tue Oct 07 11:00 AM 2008 CEST
 */

#ifndef __DSDP_WRAPPER_HPP__
#define __DSDP_WRAPPER_HPP__
#include <boost/shared_ptr.hpp>
#include <sdp_wrapper.hpp>
class DSDPWrapper: public SDPWrapper
{
	private:
		struct Impl;
		boost::shared_ptr<Impl> mImpl;

	public:
		DSDPWrapper();
		virtual ~DSDPWrapper();
		virtual AnswerT operator()(const SDPProb&);
};
#endif

