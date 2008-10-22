/*       Created   :  10/06/2008 12:54:34 AM
 *       Last Change: Wed Oct 22 11:00 PM 2008 CEST
 */

#ifndef __SDP_WRAPPER_HPP__
#define __SDP_WRAPPER_HPP__

#include <boost/numeric/ublas/fwd.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
class SDPProb;

class SDPWrapper{
	public:
		typedef boost::numeric::ublas::matrix<double> AnswerT;
		virtual AnswerT operator()(const SDPProb&);
		virtual void configure();
		virtual ~SDPWrapper();
};
#endif /* __SDP_WRAPPER_HPP__ */
