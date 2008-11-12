/*       Created   :  10/05/2008 10:25:05 PM
 *       Last Change: Sun Nov 02 09:00 PM 2008 CET
 */
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <sdp_prob.hpp>
#include <sdp_seriation_gen.hpp>
#include <sdp_seriation_prob_gen.hpp>
#include <sdp_wrapper.hpp>
#include <factory/factory.h>
#include <stats.hpp>
#include <configuration.hpp>
#include <boost/numeric/bindings/lapack/lapack.hpp>
#include <cholesky.hpp> // from third_party
#include <arg_max.hpp>
#include <boost/lambda/lambda.hpp>
#include <cmath>
#include <nana.h>
#undef C

namespace l = boost::lambda;
namespace lapack = boost::numeric::bindings::lapack;
using namespace boost::numeric;
using namespace std;

double my_sqrt(double d){return sqrt(d);}

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
	L("Generating Seriation using SDP.\n");
	const AdjMatT& adj = *adj_ptr;
	unsigned int n = adj.size1();
	I(adj.size1() == adj.size2());
	I(mSDPWrapper.get() != NULL);

	L("Generate the SDP-Problem\n");
	SDPProb prob;
	SDPSeriationProbGen sdp_probgen(adj_ptr);
	sdp_probgen(prob);

	L("Solve the SDP-Problem using an SDP Wrapper\n");
	SDPWrapper::AnswerT _X(n,n);
	noalias(_X) = (*mSDPWrapper)(prob);
	//ublas::symmetric_adaptor<SDPWrapper::AnswerT, ublas::lower> X(_X);
	SDPWrapper::AnswerT& X=_X;

//	cout << "Y = [ ";
//	cout.precision(20);
//	for(unsigned int i=0;i<X.size1();i++){
//		for(unsigned int j=0;j<X.size2();j++){
//			cout << X(i,j) << " ";
//		}
//		cout << " ; ";
//	}
//	cout << "]; \n";

#ifndef NDEBUG
	I(ublas::is_symmetric(X));
	L("Make sure we got the right thing: tr(EX)=1\n");
	using ublas::range;
	using ublas::prod;
	for(unsigned int i=0;i<prob.F.size();i++){
		ublas::matrix<double> R = prod(prob.F[i],X);
		ublas::matrix_vector_range<SDPWrapper::AnswerT> diag(R, range (0,n), range (0,n));
		double trace = ublas::sum(diag);
		//L("Trace should be: %2.3f, Trace is: %2.3f\n",prob.b[i],trace);
		I(fabs(trace-prob.b[i])<0.001);
	}
#endif

// TODO: Determine Method for rounding
#define Y_METHOD 2
#if (Y_METHOD == 1 || Y_METHOD == 3)
	// cholesky-decomposition
	L("Decompose X = V*V'\n");
	SDPWrapper::AnswerT V (X);
	if(!ulapack::chol_checked_inplace(V, true))
		throw runtime_error("Could not cholesky-decompose matrix, seems not to be positive-semidefinite.");
#endif
#if (Y_METHOD == 2 || Y_METHOD == 3)
	// eigen-decomposition
	L("Decompose X = Q * D * Q'\n");
	ublas::matrix<double,ublas::column_major> Eigv(X);
	ublas::vector<double> lambda(n);
	lapack::syev( 'V', 'L', Eigv, lambda, lapack::minimal_workspace() );
	ublas::vector<double>::iterator max_lambda = max_element(lambda.begin(),lambda.end());
	int max_lambda_idx = std::distance(lambda.begin(),max_lambda);
	ublas::vector<double> lambda_sqrt(n);
	std::transform(lambda.begin(),lambda.end(),lambda_sqrt.begin(),my_sqrt);
#endif

#if Y_METHOD == 1
	V = prod(trans(V),Eigv);
#endif

    // What do we want to maximize?
	// maximize symmetric version of C
	ublas::symmetric_adaptor<SDPProb::MatT, ublas::upper> B(prob.C);
	// maximize original C
	//SDPProb::MatT& B = prob.C;

	L("Rank reduction using random hyperplanes\n");
	ublas::vector<double> best_y(n);
	double best_y_val = -1E6;
	//srand48(gCfg().getFloat("rand_adj_mat_gen.seed"));
	ublas::vector<double> r(n);
	ublas::vector<double> y(n);
	int tries = gCfg().getInt("ser-gen.sdp-rand-plane-tries");
	ExactDescriptiveStatistics yvalstats("y_val ");
	for(int iter=0;iter<tries;iter++){
		// generate unit-vector
		generate(r.begin(), r.end(), drand48);
		r -= ublas::scalar_vector<double>(n,0.5);
		//r /= ublas::norm_2(r);

#if Y_METHOD == 1
		//use combination of cholesky-vectors
		noalias(y) = prod(ublas::trans(V),r);
		//noalias(y) = ublas::row(V,0);
		//sort(y.begin(),y.end());
		//noalias(y) = r;
#elif Y_METHOD == 2
		// use combination of eigen-vectors
		//noalias(y) = prod(Eigv,r);
		noalias(y) = sqrt(lambda(max_lambda_idx)) * ublas::column(Eigv,max_lambda_idx);
#elif Y_METHOD == 3
		// use technique from Nemirovski, Roos, Terlaky
		for(unsigned int i=0;i<r.size();i++)
			r(i) = r(i) > 0 ? 1 : -1;
		//noalias(y) =  element_prod(lambda_sqrt,prod(V,r));
		noalias(y) =  prod(V,r);
#endif
		y /= ublas::norm_2(y);

		double y_val = inner_prod(y,prod(B,y));
		yvalstats += y_val;

		if(y_val> best_y_val){
		    // the solver always _maximizes_
			best_y_val = y_val;
			best_y     = y;
		}
	}
	L("best_y_val = %2.10f\n",best_y_val);
	cout << yvalstats<<endl;

	ublas::vector<double> x = prod(sdp_probgen.mOmega_m_1_2,best_y);

//#define BEST_ELEM(X,tmp) util::arg_max(X.begin(),X.end(),tmp, l::_1*l::_1)
#define BEST_ELEM(X,tmp) max_element(X.begin(),X.end())

	// tricky: make sure x > 0 at all times.
	x += ublas::scalar_vector<double>(x.size(), *min_element(x.begin(),x.end()) + 1);

	SeriationT ret;
	std::vector<bool> done(x.size(),false);

	// find highest component of x
	double tmp_d;
	ublas::vector<double>::iterator it = BEST_ELEM(x,tmp_d);
	int idx = std::distance(x.begin(),it);

	L("Determine Actual Path through Graph.\n");
	for(unsigned int i=0;i<n;i++){
		// mark as visited
		done[idx] = false;

		// make sure we do not visit again
		*it = 0;

		// [m,i] = max(x.*A(:,i));
		ublas::vector<double> tmp = ublas::element_prod(x,ublas::column(adj,idx));
		it = BEST_ELEM(tmp,tmp_d);

		int old_idx = idx;
		idx = std::distance(tmp.begin(),it);

		// point it in x, not tmp:
		it = x.begin();
		it += idx;

		if( (done[idx] || adj(old_idx,idx)<0.5) && ret.size()<n)
		{
			// we reached a dead end
			double tmp;
			it = BEST_ELEM(x,tmp_d);
			idx = std::distance(x.begin(),it);
			//L("jumping from %d to %d\n", old_idx, idx );
			//break;
		}else{
			// remember transition
			ret.push_back(make_pair(old_idx,idx));
		}
	}

	//for(SeriationT::iterator si=ret.begin();si!=ret.end();si++){
		//cout << "path("<<si->first<<", "<<si->second<<")."<<endl;
	//}

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

