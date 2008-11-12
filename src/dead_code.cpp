/****************************************************
 * Taken from main after realizing hancock is bad
 ****************************************************/

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
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
