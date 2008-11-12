/*       Created   :  10/05/2008 10:21:12 PM
 *       Last Change: Wed Nov 12 07:00 PM 2008 CET
 */

#ifndef __SDP_SERIATION_GEN_HPP__
#define __SDP_SERIATION_GEN_HPP__
#include <memory>
#include <seriation_gen.hpp>

class SDPWrapper;

class SDPSeriationGen : public SeriationGen
{
	private:
		struct Impl;
		boost::shared_ptr<Impl> mImpl;
	public:
        SDPSeriationGen();
		void setSDPWrapper(std::auto_ptr<SDPWrapper> w);
		virtual SeriationT operator()(boost::shared_ptr<AdjMatT>);
		virtual void configure();
		virtual ~SDPSeriationGen();
};
#endif /* __SDP_SERIATION_GEN_HPP__ */

