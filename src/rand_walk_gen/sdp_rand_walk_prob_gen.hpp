/*       Created   :  10/06/2008 12:23:53 AM
 *       Last Change: Tue Oct 07 12:00 AM 2008 CEST
 */

#ifndef __SDP_RAND_WALK_PROB_GEN_HPP__
#define __SDP_RAND_WALK_PROB_GEN_HPP__
#include <sdp_prob.hpp>
class SDPRandWalkProbGen{
	public:
		typedef boost::numeric::ublas::matrix<double,
				boost::numeric::ublas::column_major> AdjMatT;
		const AdjMatT&  mAdjMat;
		SDPRandWalkProbGen(const AdjMatT&);
		void operator()(SDPProb&);
};
#endif /* __SDP_RAND_WALK_PROB_GEN_HPP__ */

