/*       Created   :  10/06/2008 12:57:46 AM
 *       Last Change: Mon Oct 06 01:00 AM 2008 CEST
 */

#ifndef __CSDP_WRAPPER_HPP__
#define __CSDP_WRAPPER_HPP__
#include <sdp_wrapper.hpp>
class SDPAWrapper: public SDPWrapper
{
    public:
        virtual AnswerT operator()(const SDPProb&);
        virtual ~SDPAWrapper();
	private:
        void writeSDPAInputFile(const SDPProb&, const char*);
        void runSDPA(const char* in, const char* out,const char* par);
};
#endif

