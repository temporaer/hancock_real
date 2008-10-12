/*       Created   :  10/06/2008 12:54:34 AM
 *       Last Change: Tue Oct 07 08:00 PM 2008 CEST
 */

#ifndef __SDP_WRAPPER_HPP__
#define __SDP_WRAPPER_HPP__

#include <boost/numeric/ublas/fwd.hpp>
//#include <sdp_prob.hpp>
class SDPProb;

class SDPWrapper{
	public:
		typedef boost::numeric::ublas::vector<double> AnswerT;
		virtual AnswerT operator()(const SDPProb&);
		virtual ~SDPWrapper();
};
#endif /* __SDP_WRAPPER_HPP__ */
