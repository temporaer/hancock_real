/*       Created   :  10/05/2008 10:25:05 PM
 *       Last Change: Tue Oct 07 09:00 PM 2008 CEST
 */
#include <nana.h>
#include <sdp_seriation_gen.hpp>
#include <sdp_seriation_prob_gen.hpp>
#include <sdp_wrapper.hpp>
#include <factory/factory.h>

using namespace boost::numeric::ublas;
using namespace std;
SeriationGen::SeriationT SDPSeriationGen::operator()(const AdjMatT& adj)
{
	I(adj.size1() == adj.size2());

	// Generate the SDP-Problem
	SDPProb prob;
	SDPSeriationProbGen sdp_probgen(adj);
	sdp_probgen(prob);

	// Solve the SDP-Problem using an SDP Wrapper
	(*mSDPWrapper)(prob);

	SeriationT ret;
	return ret;
}


