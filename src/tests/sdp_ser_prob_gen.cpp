#define BOOST_TEST_MODULE sdp_ser_gen_test

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <sdp_prob.hpp>
#include <sdp_seriation_prob_gen.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/shared_ptr.hpp>
using namespace std;
using namespace boost;

struct Fixture{
	const int n;
	shared_ptr<SDPSeriationProbGen> sdpserprobgen;
	shared_ptr<SDPSeriationProbGen::AdjMatT> am;
	
	Fixture()
		:  n(6)
		 , am(new SDPSeriationProbGen::AdjMatT(n,n))
	{
		for(int i=0;i<n;i++)
			for(int j=0;j<n;j++)
			{
				(*am)(i,j) = ((i+j)%2>0)?1.0:0.0;
			}
		sdpserprobgen.reset(new SDPSeriationProbGen(am));
	}
	~Fixture(){
	}
};

BOOST_FIXTURE_TEST_SUITE( suite, Fixture )


BOOST_AUTO_TEST_CASE( testF1 )
{
	SDPProb prob;
	(*sdpserprobgen)(prob);
	BOOST_REQUIRE_GE(prob.F.size(),1);
	BOOST_REQUIRE_EQUAL(prob.F[0].size1(),n);
	BOOST_REQUIRE_EQUAL(prob.F[0].size2(),n);
	boost::numeric::ublas::identity_matrix<double> im(n);
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++)
			BOOST_CHECK_EQUAL(im(i,j),prob.F[0](i,j));

}

BOOST_AUTO_TEST_CASE( testC )
{
	double d[6][6] = 
	{
		// calculated with octave, omega^(1/2)*A*omega^(-1/2)
		{0.00000,  0.70711,  0.00000,  0.70711,  0.00000,  1.00000},
		{1.41421,  0.00000,  1.00000,  0.00000,  1.00000,  0.00000},
		{0.00000,  1.00000,  0.00000,  1.00000,  0.00000,  1.41421},
		{1.41421,  0.00000,  1.00000,  0.00000,  1.00000,  0.00000},
		{0.00000,  1.00000,  0.00000,  1.00000,  0.00000,  1.41421},
		{1.00000,  0.00000,  0.70711,  0.00000,  0.70711,  0.00000},
	};
	SDPProb prob;
	(*sdpserprobgen)(prob);
	BOOST_REQUIRE_EQUAL(prob.C.size1(),n);
	BOOST_REQUIRE_EQUAL(prob.C.size2(),n);
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++)
			BOOST_CHECK_CLOSE(d[i][j],prob.C(i,j),0.01);

}

BOOST_AUTO_TEST_CASE( testOmegas )
{
	SDPProb prob;
	(*sdpserprobgen)(prob);
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++)
		{
			if(i!=j){
				BOOST_REQUIRE_EQUAL(sdpserprobgen->mOmega_1_2(i,j),0);
				BOOST_REQUIRE_EQUAL(sdpserprobgen->mOmega_m_1_2(i,j),0);
				continue;
			}
			if(i==n-1 || i == 0){
				BOOST_REQUIRE_EQUAL(sdpserprobgen->mOmega_1_2(i,j),1);
				BOOST_REQUIRE_EQUAL(sdpserprobgen->mOmega_m_1_2(i,j),1);
				continue;
			}
			BOOST_REQUIRE_EQUAL(sdpserprobgen->mOmega_1_2(i,j),sqrt(2));
			BOOST_REQUIRE_EQUAL(sdpserprobgen->mOmega_m_1_2(i,j),1/sqrt(2));

		}
	
}

BOOST_AUTO_TEST_CASE( test_b )
{
	SDPProb prob;
	(*sdpserprobgen)(prob);
	BOOST_REQUIRE_EQUAL(prob.b.size(),1);
	BOOST_CHECK_EQUAL(prob.b[0],1);
}


BOOST_AUTO_TEST_SUITE_END()

