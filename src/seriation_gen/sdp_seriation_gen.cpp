/*       Created   :  10/05/2008 10:25:05 PM
 *       Last Change: Thu Oct 23 10:00 AM 2008 CEST
 */
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <sdp_prob.hpp>
#include <sdp_seriation_gen.hpp>
#include <sdp_seriation_prob_gen.hpp>
#include <sdp_wrapper.hpp>
#include <factory/factory.h>
#include <nana.h>
#include <cholesky.hpp> // from third_party

using namespace boost::numeric;
using namespace std;

// Implementation of SeriationGen

struct SDPSeriationGen::Impl{
    std::auto_ptr<SDPWrapper> mSDPWrapper;
    typedef SeriationGen::AdjMatT AdjMatT;
    SeriationGen::SeriationT operator()(boost::shared_ptr<AdjMatT> adj);
    void setSDPWrapper(std::auto_ptr<SDPWrapper> w);
    ~Impl();
};
SDPSeriationGen::Impl::~Impl(){
}

SeriationGen::SeriationT SDPSeriationGen::Impl::operator()(boost::shared_ptr<AdjMatT> adj_ptr)
{
	AdjMatT& adj = *adj_ptr;
	I(adj.size1() == adj.size2());
	I(mSDPWrapper.get() != NULL);

	// Generate the SDP-Problem
	SDPProb prob;
	SDPSeriationProbGen sdp_probgen(adj_ptr);
	sdp_probgen(prob);

	// Solve the SDP-Problem using an SDP Wrapper
	SDPWrapper::AnswerT X = (*mSDPWrapper)(prob);
	
	// Make sure we got the right thing: tr(EX)=1
#ifndef NDEBUG
	using ublas::range;
	using ublas::prod;
	int n = X.size1();
	for(unsigned int i=0;i<prob.F.size();i++){
		ublas::matrix<double> R = prod(prob.F[i],X);
		ublas::matrix_vector_range<SDPWrapper::AnswerT> diag(R, range (0,n), range (0,n));
		double trace = ublas::sum(diag); 
		L("Trace should be: %2.3f, Trace is: %2.3f\n",prob.b[i],trace);
		I(fabs(trace-prob.b[i])<0.001);
	}
#endif

	// decompose X = V*V'
	SDPWrapper::AnswerT& V = X;
	ulapack::chol_checked_inplace(V, true);


	// random hyperplane technique

	SeriationT ret;
	return ret;
}


void SDPSeriationGen::Impl::setSDPWrapper(std::auto_ptr<SDPWrapper> w){
  mSDPWrapper = w;
}




// Wrap SeriationGen::Impl

SDPSeriationGen::SDPSeriationGen()
:mImpl(new Impl)
{

}


SeriationGen::SeriationT SDPSeriationGen::operator()(boost::shared_ptr<AdjMatT> adj)
{
    return (*mImpl)(adj);
}

void SDPSeriationGen::setSDPWrapper(std::auto_ptr<SDPWrapper> w){
    mImpl->setSDPWrapper(w);
}

SDPSeriationGen::~SDPSeriationGen()
{
L("Destroying SDPSeriationGen\n");
}

