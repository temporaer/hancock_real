/*       Created   :  10/03/2008 08:22:01 PM
 *       Last Change: Mon Oct 06 06:00 PM 2008 CEST
 */
#include <string>
#include <iostream>
#include <nana.h>
#include "configuration.hpp"
#include <sdp_rand_walk_gen.hpp>

#include <boost/numeric/ublas/matrix.hpp>
using namespace boost::numeric::ublas;
using namespace std;

Configuration gCfg;

using namespace std;

int main(int argc, char* argv[]){
	gCfg.parsecfg(argc,argv);

	typedef matrix<double,row_major>    AdjMatRT;
	typedef matrix<double,column_major> AdjMatCT;

	int n = 4;
	AdjMatRT adjmat(n,n);
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++)
			adjmat(i,j) = i*10+j;

	AdjMatCT adjmat2(adjmat);

	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++)
			I(adjmat(i,j) == adjmat2(i,j));

	double* p = (double*)(&adjmat2(0,0));
	for(int i=0;i<n*n;i++)
		cout<<p[i]<<endl;


	SDPRandWalkGen walkgen;

	RandWalkGen::RandWalkT randwalk  = walkgen(adjmat);
}
