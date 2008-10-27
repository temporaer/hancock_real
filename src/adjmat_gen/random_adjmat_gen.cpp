#include <random_adjmat_gen.hpp>
#include <iostream>
#include <boost/numeric/ublas/matrix.hpp>
#include <configuration.hpp>

using namespace boost;
using namespace std;
void RandomAdjMatGen::configure()
{
	this->setMatrixSize(gCfg().getInt("rand_adj_mat_gen.size"));
	this->setConnectionProb(gCfg().getFloat("rand_adj_mat_gen.prob"));
	this->setWeightedEdges(gCfg().getBool("rand_adj_mat_gen.weighted"));
	this->setSeed(gCfg().getFloat("rand_adj_mat_gen.seed"));
}

shared_ptr<RandomAdjMatGen::AdjMatT> RandomAdjMatGen::operator()()
{
	shared_ptr<AdjMatT> adjmat_ptr(new AdjMatT(mMatrixSize,mMatrixSize));
	
	AdjMatT& adjmat = *adjmat_ptr;
	int n = adjmat.size1();
	srand48(mSeed);
	for(int i=0;i<n;i++)
        for(int j=i;j<n;j++){
			if(mWeightedEdges){
				if(drand48()<mConnectionProb)
					adjmat(i,j) = 0;
				else
					adjmat(i,j) = drand48();
			}else{
				adjmat(i,j) = drand48()>mConnectionProb;
			}
            adjmat(j,i) = adjmat(i,j);
		}
	
	cout<<"A = [ ";
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			cout << adjmat(i,j) <<" ";
		}
		cout << "; ";
	}
	cout << " ]; \n";
	return adjmat_ptr;
}

