/*       Created   :  10/06/2008 12:36:07 AM
 *       Last Change: Thu Oct 23 10:00 AM 2008 CEST
 */

#include <sdp_seriation_prob_gen.hpp>
#include <sdp_prob.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <nana.h>
#undef C

namespace ublas = boost::numeric::ublas;
using namespace std;

SDPSeriationProbGen::SDPSeriationProbGen(boost::shared_ptr<AdjMatT>a)
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
	
	Matrix& adj = *mAdjMat;

	I(adj.size1() == adj.size2());

	int n = adj.size1();

	calcOmega(n);

	// C = F_0
	prob.C = prod(mOmega_1_2, Matrix(prod(adj,mOmega_m_1_2)));

#ifndef NDEBUG
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++){
			I(mOmega_1_2(i,j) == mOmega_1_2(i,j));
			I(mOmega_m_1_2(i,j) == mOmega_m_1_2(i,j));
			I(adj(i,j) == adj(i,j));
			I(prob.C(i,j) == prob.C(i,j));
		}
#endif

	// F_1
	identity_matrix<double> Id(n);
	prob.F.push_back(Id);

	// b_1
	prob.b = ublas::vector<double>(1);
	prob.b(0) = 1;
}

