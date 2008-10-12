/*       Created   :  10/06/2008 12:23:53 AM
 *       Last Change: Tue Oct 07 09:00 PM 2008 CEST
 */

#ifndef __SDP_SERIATION_PROB_GEN_HPP__
#define __SDP_SERIATION_PROB_GEN_HPP__
#include <sdp_prob.hpp>
class SDPSeriationProbGen{
	public:
		typedef boost::numeric::ublas::matrix<double,
				boost::numeric::ublas::column_major> AdjMatT;
		SDPSeriationProbGen(const AdjMatT&);
		void operator()(SDPProb&);

    private:
        const AdjMatT&  mAdjMat;
		AdjMatT         mOmega;
		AdjMatT         mOmega_1_2;
		AdjMatT         mOmega_m_1_2;
		void calcOmega(int n);
};
#endif /* __SDP_SERIATION_PROB_GEN_HPP__ */

