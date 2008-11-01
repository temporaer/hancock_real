/*       Created   :  10/03/2008 08:22:01 PM
 *       Last Change: Fri Oct 31 05:00 PM 2008 CET
 */

#include <dlfcn.h>
#include <string>
#include <iostream>
#include <fstream>
#include <factory/factory.h>
#include "configuration.hpp"

#include <sdp_seriation_gen.hpp>
#include <sdp_wrapper.hpp>

#include <random_adjmat_gen.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/shared_ptr.hpp>
#include <nana.h>
#undef C
using namespace boost::numeric::ublas;
using namespace boost;
using namespace std;

std::map<int, int> idxmap;

template<class T, class M>
ofstream& pathForPerl(ofstream& o,
		const M& adj,
		const T& path,
		bool is_perm){
	for(int i=0;i<path.size();i++){
		if(!is_perm){
			o << path[i].first << " " << path[i].second << " ";
		}
		else{
			I(adj(idxmap[path[i].first],idxmap[path[i].second])>0.5);
			o << idxmap[path[i].first] << " " << idxmap[path[i].second] << " ";
		}
	}
	o << endl;
	return o;
}
template<class T,class U>
void writeToFile(
		const char* fn,
		const T&    adj,
		const U&    path,
		const U&    path_perm
		){
	ofstream o(fn);
	o << "A = [ " ;
	for(int row=0;row<adj.size2();row++){
		for(int col=0;col<adj.size1();col++){
			o << adj(row,col) << " ";
		}
		if(row != adj.size2()-1)
			o << " ; ";
	}
	o << " ];"<<endl;

	o << "path = [ ";
	for(int i=0;i<path.size();i++)
	{
		o << path[i].first << " ";
	}
	if(path.size()>0)
		o << path.back().second ;
	o << " ] ; "<< endl;

	o << "path_perm = [ ";
	for(int i=0;i<path_perm.size();i++)
	{
		o << idxmap[path_perm[i].first] << " ";
	}
	if(path_perm.size()>0) 
		o << idxmap[path_perm.back().second];
	o << " ] ; "<< endl;
}

shared_ptr<RandomAdjMatGen::AdjMatT> copyAndPermute(shared_ptr<RandomAdjMatGen::AdjMatT> org_ptr){
	L("Generating random permutation...");
	RandomAdjMatGen::AdjMatT& org = *org_ptr;
	int n = org.size1();

	std::vector<int> idxs;
	for(int i=0;i<n;i++)
		idxs.push_back(i);

	int r = (int) pow(2,drand48()*n);
	std::vector<int>::iterator start=idxs.begin(),end=idxs.end();
	random_shuffle (idxs.begin(), idxs.end());

	cout << "Random Permutation: ";
	for(int i=0;i<n;i++){
		cout << idxs[i]<<" ";
		idxmap[idxs[i]] = i;
	}
	cout << endl;

	// create new matrix
	shared_ptr<RandomAdjMatGen::AdjMatT> ret_ptr(new RandomAdjMatGen::AdjMatT(org.size1(),org.size2())); 
	RandomAdjMatGen::AdjMatT& ret = *ret_ptr;


	// fill with values
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++)
			ret(i,j) = org(idxmap[i],idxmap[j]);

	return ret_ptr;
}

template<class T,class U>
void comparePaths(shared_ptr<T> adj_ptr, const U& p, const U& p_perm, bool check=true){
	T& adj = *adj_ptr;
	int nfound = 0;
	int nnotfound =0;
	for(int i=0;i<p.size();i++){
		bool found = false;
		for(int j=0;j<p_perm.size();j++){
			if(check){
				I(adj(idxmap[p_perm[j].first],idxmap[p_perm[j].second])>0.5);
			}
			if(
					(p[i].first  == idxmap[p_perm[j].first]  && p[i].second == idxmap[p_perm[j].second])
				||  (p[i].first  == idxmap[p_perm[j].second] && p[i].second == idxmap[p_perm[j].first])
					){
				found = true;
				break;
			}
		}
		if(found)
			nfound++;
		else
			nnotfound++;
	}
	L("Of %d/%d pairs, %d were found in p and permutation.\n", p.size(),p_perm.size(), nfound);
}

template<class T,class U>
void comparePathsPlain(shared_ptr<T> adj_ptr, const U& p, const U& p_perm){
	T& adj = *adj_ptr;
	int nfound = 0;
	int nnotfound =0;
	for(int i=0;i<p.size();i++){
		bool found = false;
		for(int j=0;j<p_perm.size();j++){
			if(
					(p[i].first  == p_perm[j].first  && p[i].second == p_perm[j].second)
				||  (p[i].first  == p_perm[j].second && p[i].second == p_perm[j].first)
			){
				found = true;
				break;
			}
		}
		if(found)
			nfound++;
		else
			nnotfound++;
	}
	L("Of %d/%d pairs, %d were found in p and repetition.\n", p.size(),p_perm.size(), nfound);
}

int main(int argc, char* argv[]){

	dlopen("sdp_wrappers/libsdp_wrappers.so",RTLD_LAZY);
	
	gCfg().parsecfg(argc,argv);

	RandomAdjMatGen matgen;
	matgen.configure();
	shared_ptr<RandomAdjMatGen::AdjMatT> adjmat_ptr(matgen());
	shared_ptr<RandomAdjMatGen::AdjMatT> perm_ptr  (copyAndPermute(adjmat_ptr));

	string sdp_wrapper_name          = gCfg().getString("ser-gen.sdp-wrapper");
	auto_ptr<SDPWrapper> sdp_wrapper = genericFactory<SDPWrapper>::instance().create(sdp_wrapper_name);
	sdp_wrapper->configure();

	SDPSeriationGen walkgen;
	walkgen.setSDPWrapper( sdp_wrapper );

	SeriationGen::SeriationT randwalk,
		randwalk_again1,
		randwalk_again2,
		randwalk_perm; 

#if 0
	try{ randwalk        = walkgen(adjmat_ptr);
	}catch(const exception& e){ L("generating randwalk failed: %s.\n",e.what()); }
	try{ randwalk_again1 = walkgen(adjmat_ptr);
	}catch(const exception& e){ L("generating randwalk_again1 failed: %s.\n",e.what()); }
	try{ randwalk_again2 = walkgen(adjmat_ptr);
	}catch(const exception& e){ L("generating randwalk_again2 failed: %s.\n",e.what()); }
	try{ randwalk_perm   = walkgen(perm_ptr);
	}catch(const exception& e){ L("generating randwalk_perm failed: %s.\n",e.what()); }


	comparePaths(adjmat_ptr, randwalk, randwalk_perm);
	comparePathsPlain(adjmat_ptr, randwalk, randwalk_again1);
	comparePathsPlain(adjmat_ptr, randwalk, randwalk_again2);
	comparePathsPlain(adjmat_ptr, randwalk_again1, randwalk_again2);

	const char* out_file = gCfg().getString("output").c_str();
	writeToFile(out_file,*adjmat_ptr,randwalk, randwalk_perm);
#else
	ofstream paths("paths.dat");
	for(int i=0;i<20;i++){
		shared_ptr<RandomAdjMatGen::AdjMatT> perm  (copyAndPermute(adjmat_ptr));
		try{ randwalk_perm   = walkgen(perm);
		}catch(const exception& e){ L("generating randwalk_perm failed: %s.\n",e.what()); }
		pathForPerl(paths,*adjmat_ptr,randwalk_perm, true);
	}
#endif
}
