/*       Created   :  10/05/2008 10:25:05 PM
 *       Last Change: Tue Oct 07 12:00 AM 2008 CEST
 */
#include <nana.h>
#include <sdp_rand_walk_gen.hpp>
#include <sdp_rand_walk_prob_gen.hpp>

using namespace boost::numeric::ublas;
using namespace std;
RandWalkGen::RandWalkT SDPRandWalkGen::operator()(const AdjMatT& adj)
{
	I(adj.size1() == adj.size2());

	// Generate the SDP-Problem
	SDPProb prob;
	SDPRandWalkProbGen sdp_probgen(adj);
	sdp_probgen(prob);

	// Solve the SDP-Problem using an SDP Wrapper
	RandWalkT ret;
	return ret;
}


