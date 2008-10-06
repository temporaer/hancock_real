/*       Created   :  10/05/2008 10:25:05 PM
 *       Last Change: Mon Oct 06 06:00 PM 2008 CEST
 */
#include <sdp_rand_walk_gen.hpp>
#include <sdp_rand_walk_prob_gen.hpp>
using namespace boost::numeric::ublas;
using namespace std;
RandWalkGen::RandWalkT SDPRandWalkGen::operator()(const AdjMatT& adj)
{
	matrix<double> omega;

	SDPRandWalkProbGen sdp_probgen(omega,adj);

	SDPProb prob = sdp_probgen();

	RandWalkT ret;
	return ret;
}


