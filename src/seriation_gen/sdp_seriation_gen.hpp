/*       Created   :  10/05/2008 10:21:12 PM
 *       Last Change: Tue Oct 07 09:00 PM 2008 CEST
 */

#ifndef __SDP_SERIATION_GEN_HPP__
#define __SDP_SERIATION_GEN_HPP__
#include <memory>
#include <seriation_gen.hpp>
#include <sdp_wrapper.hpp>
class SDPSeriationGen : public SeriationGen 
{
	private:
		std::auto_ptr<SDPWrapper> mSDPWrapper;
	public:
		inline  void setSDPWrapper(std::auto_ptr<SDPWrapper> w){mSDPWrapper = w;}
		virtual SeriationT operator()(const AdjMatT&);
};
#endif /* __SDP_SERIATION_GEN_HPP__ */

