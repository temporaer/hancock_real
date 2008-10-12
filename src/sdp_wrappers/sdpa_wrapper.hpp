/*       Created   :  10/06/2008 12:57:46 AM
 *       Last Change: Sun Oct 12 10:00 PM 2008 CEST
 */

#ifndef __CSDP_WRAPPER_HPP__
#define __CSDP_WRAPPER_HPP__
#include <sdp_wrapper.hpp>
class SDPAWrapper: public SDPWrapper
{
    public:
        //! solve an @see SDPProb SDP-problem
        virtual AnswerT operator()(const SDPProb&);
        virtual ~SDPAWrapper();
	private:

        /*! Write an input file for sdpa to read.
         *  \param p   the problem to print
         *  \param in  the file to write to */
        void writeSDPAInputFile(const SDPProb& p, const char* in);

        /*! Read the output of sdpa.
         *  \param out the name of the sdpa-output file to read
         *  \param ret the vector where we should record the answer */
        bool readSDPAOutputFile(const char* out, AnswerT& ret);

        /*! Run the sdpa-binary.
         *  \param in where the problem was written to
         *  \param out where the solution shall be written to
         *  \param par the parameter file for sdpa (see its doc) */
        void runSDPA(const char* in, const char* out,const char* par);
};
#endif

