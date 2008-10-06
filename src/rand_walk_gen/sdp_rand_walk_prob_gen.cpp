/*       Created   :  10/06/2008 12:36:07 AM
 *       Last Change: Mon Oct 06 12:00 AM 2008 CEST
 */

#include <sdp_rand_walk_prob_gen.hpp>

SDPRandWalkProbGen::SDPRandWalkProbGen(const OmegaT&o, const AdjMatT&a)
	: mOmega(o)
	 ,mAdjMat(a)
{
}


SDPProb SDPRandWalkProbGen::operator()()
{
	SDPProb prob;
	return prob;
}

