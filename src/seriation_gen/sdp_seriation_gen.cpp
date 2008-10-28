/*       Created   :  10/05/2008 10:25:05 PM
 *       Last Change: Tue Oct 28 01:00 AM 2008 CET
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

//	cout << "Y = [ ";
//	for(unsigned int i=0;i<X.size1();i++){
//		for(unsigned int j=0;j<X.size2();j++){
//			cout << X(i,j) << " ";
//		}
//		cout << " ; ";
//	}
//	cout << "]; \n";
	
	unsigned int n = X.size1();
	// Make sure we got the right thing: tr(EX)=1
#ifndef NDEBUG
	using ublas::range;
	using ublas::prod;
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
	if(!ulapack::chol_checked_inplace(V, true))
		throw runtime_error("Could not cholesky-decompose matrix, seems not to be positive-semidefinite.");

//	cout << "V = [ ";
//	for(unsigned int i=0;i<V.size1();i++){
//		for(unsigned int j=0;j<V.size2();j++){
//			cout << V(i,j) << " ";
//		}
//		cout << " ; ";
//	}
//	cout << "]; \n";

	// random hyperplane technique
	ublas::vector<double> min_y;
	float min_y_val = 1E6;
	// TODO: intelligenter abbrechen
	for(int iter=0;iter<10000;iter++){

		// generate unit-vector
		ublas::vector<double> r(n);
		generate(r.begin(), r.end(), rand); 
		r /= ublas::norm_2(r);

		ublas::vector<double> y = prod(ublas::trans(V),r);
		sort(y.begin(),y.end());

		float y_val = inner_prod(r,prod(prob.C,r));
		if(y_val< min_y_val){
			min_y_val = y_val;
			min_y     = y;
		}
	}

	ublas::vector<double> x = prod(sdp_probgen.mOmega_m_1_2,min_y);

	// tricky: make sure x > 0 at all times.
	x += ublas::scalar_vector<double>(x.size(), *min_element(x.begin(),x.end()) + 1);

	SeriationT ret;
	std::vector<bool> done(x.size(),false);

	// find highest component of x
	ublas::vector<double>::iterator it = max_element(x.begin(),x.end());
	int idx = std::distance(x.begin(),it);

	for(unsigned int i=0;i<n;i++){
		// mark as visited
		done[idx] = false;

		// make sure we do not visit again
		*it = -1;

		// [m,i] = max(x.*A(:,i));
		ublas::vector<double> tmp = ublas::element_prod(x,ublas::column(adj,idx));
		it = max_element(tmp.begin(),tmp.end());

		int old_idx = idx;
		idx = std::distance(tmp.begin(),it);

		// point it in x, not tmp:
		it = x.begin();
		it += idx;

		if( (done[idx] || adj(old_idx,idx)<0.5) && ret.size()<n)
		{
			// we reached a dead end
			it  = max_element(x.begin(),x.end());
			idx = std::distance(x.begin(),it);
			L("jumping from %d to %d\n", old_idx, idx );
		}else{
			// remember transition
			ret.push_back(make_pair(old_idx,idx));
		}
	}

	for(SeriationT::iterator si=ret.begin();si!=ret.end();si++){
		cout << "path("<<si->first<<", "<<si->second<<")."<<endl;
	}

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
}

