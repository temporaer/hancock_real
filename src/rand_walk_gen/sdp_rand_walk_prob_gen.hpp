/*       Created   :  10/06/2008 12:23:53 AM
 *       Last Change: Mon Oct 06 12:00 AM 2008 CEST
 */

#ifndef __SDP_RAND_WALK_PROB_GEN_HPP__
#define __SDP_RAND_WALK_PROB_GEN_HPP__
#include <sdp_prob.hpp>
class SDPRandWalkProbGen{
	public:
		typedef boost::numeric::ublas::matrix<double> AdjMatT;
		typedef boost::numeric::ublas::matrix<double> OmegaT;
		const OmegaT&   mOmega;
		const AdjMatT&  mAdjMat;
		SDPRandWalkProbGen(const OmegaT&, const AdjMatT&);
		SDPProb operator()();
};
#endif /* __SDP_RAND_WALK_PROB_GEN_HPP__ */

