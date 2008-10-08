/*       Created   :  10/05/2008 10:21:12 PM
 *       Last Change: Tue Oct 07 09:00 PM 2008 CEST
 */

#ifndef __SDP_SERIATION_GEN_HPP__
#define __SDP_SERIATION_GEN_HPP__
#include <memory>
#include <boost/shared_ptr.hpp>
#include <seriation_gen.hpp>
//#include <sdp_wrapper.hpp>

class SDPWrapper;

class SDPSeriationGen : public SeriationGen
{
	private:
		struct Impl;
		boost::shared_ptr<Impl> mImpl;
	public:
        SDPSeriationGen();
		void setSDPWrapper(std::auto_ptr<SDPWrapper> w);
		virtual SeriationT operator()(const AdjMatT&);
		virtual ~SDPSeriationGen();
};
#endif /* __SDP_SERIATION_GEN_HPP__ */

