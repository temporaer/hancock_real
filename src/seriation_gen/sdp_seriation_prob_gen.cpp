/*       Created   :  10/06/2008 12:36:07 AM
 *       Last Change: Thu Oct 16 04:00 PM 2008 CEST
 */

#include <nana.h>
#include <sdp_seriation_prob_gen.hpp>
#include <sdp_prob.hpp>
//#include <boost/numeric/ublas/io.hpp>

namespace ublas = boost::numeric::ublas;
using namespace std;

SDPSeriationProbGen::SDPSeriationProbGen(const AdjMatT&a)
	: mAdjMat(a)
{
}

void SDPSeriationProbGen::calcOmega(int n){
    I(n>0);
    using ublas::identity_matrix;
	typedef SDPProb::MatT Matrix ;

//	// omega
//	Matrix mOmega( 2*identity_matrix<double>(n) );
//	mOmega(0,0)     = 1;
//	mOmega(n-1,n-1) = 1;

	// omega^(1/2)
	mOmega_1_2 = Matrix( sqrt(2) * identity_matrix<double>(n) );
	mOmega_1_2(0,0)     = 1;
	mOmega_1_2(n-1,n-1) = 1;

	// omega^(-1/2)
	mOmega_m_1_2 = Matrix( (1.0/sqrt(2)) * identity_matrix<double>(n) );
	mOmega_m_1_2(0,0)     = 1;
	mOmega_m_1_2(n-1,n-1) = 1;
}

void SDPSeriationProbGen::operator()(SDPProb& prob)
{
	using ublas::identity_matrix;
	typedef SDPProb::MatT Matrix ;

	I(mAdjMat.size1() == mAdjMat.size2());

	int n = mAdjMat.size1();

	calcOmega(n);

	// C = F_0
	prob.C = prod(mOmega_1_2, Matrix(prod(mAdjMat,mOmega_m_1_2)));

	// F_1
	identity_matrix<double> Id(n);
	prob.F.push_back(Id);

	// b_1
	prob.b = ublas::vector<double>(1);
	prob.b(0) = 1;
}

