/*       Created   :  10/06/2008 12:23:53 AM
 *       Last Change: Tue Oct 28 01:00 AM 2008 CET
 */

#ifndef __SDP_SERIATION_PROB_GEN_HPP__
#define __SDP_SERIATION_PROB_GEN_HPP__
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/shared_ptr.hpp>
class SDPProb;
class SDPSeriationProbGen{
	public:
		typedef boost::numeric::ublas::diagonal_matrix<double> OmegaT;
		typedef boost::numeric::ublas::matrix<double,
				boost::numeric::ublas::column_major> AdjMatT;
		SDPSeriationProbGen(boost::shared_ptr<AdjMatT>);
		void operator()(SDPProb&);

    public:
		boost::shared_ptr<AdjMatT>  mAdjMat;
		OmegaT         mOmega;
		OmegaT         mOmega_1_2;
		OmegaT         mOmega_m_1_2;
		void calcOmega(int n);
};
#endif /* __SDP_SERIATION_PROB_GEN_HPP__ */

