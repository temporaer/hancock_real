/*       Created   :  10/06/2008 12:36:07 AM
 *       Last Change: Tue Oct 07 12:00 AM 2008 CEST
 */

#include <nana.h>
#include <sdp_rand_walk_prob_gen.hpp>

namespace ublas = boost::numeric::ublas;
using namespace std;

SDPRandWalkProbGen::SDPRandWalkProbGen(const AdjMatT&a)
	: mAdjMat(a)
{
}


void SDPRandWalkProbGen::operator()(SDPProb& prob)
{
	using ublas::identity_matrix;
	typedef AdjMatT Matrix ;

	I(mAdjMat.size1() == mAdjMat.size2());

	int n = mAdjMat.size1();

	// omega
	Matrix omega( 2*identity_matrix<double>(n) );
	omega(0,0)     = 1;
	omega(n-1,n-1) = 1;

	// omega^(1/2)
	Matrix omega_1_2 ( sqrt(2) * identity_matrix<double>(n) );
	omega_1_2(0,0)     = 1;
	omega_1_2(n-1,n-1) = 1;

	// omega^(-1/2)
	Matrix omega_m1_2 ( (1.0/sqrt(2)) * identity_matrix<double>(n) );
	omega_m1_2(0,0)     = 1;
	omega_m1_2(n-1,n-1) = 1;

	// C
	noalias(prob.C) = prod(omega_1_2, Matrix(prod(mAdjMat,omega_m1_2)));

	// F_0
	identity_matrix<double> Id(n);
	prob.F.push_back(Id);

	// b_0
	prob.b = ublas::vector<double>(1);
	prob.b(0) = 1;
}

