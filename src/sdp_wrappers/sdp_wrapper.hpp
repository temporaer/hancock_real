/*       Created   :  10/06/2008 12:54:34 AM
 *       Last Change: Mon Oct 06 01:00 AM 2008 CEST
 */

#ifndef __SDP_WRAPPER_HPP__
#define __SDP_WRAPPER_HPP__

#include <boost/numeric/ublas/matrix.hpp>
#include <sdp_prob.hpp>

class SDPWrapper{
	public:
		typedef boost::numeric::ublas::vector<double> AnswerT;
		virtual AnswerT operator()(const SDPProb&)=0;
};
#endif /* __SDP_WRAPPER_HPP__ */
