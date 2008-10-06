/*       Created   :  10/05/2008 07:45:51 PM
 *       Last Change: Sun Oct 05 10:00 PM 2008 CEST
 */

#ifndef __RAND_WALK_GEN_HPP__
#define __RAND_WALK_GEN_HPP__
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>

class RandWalkGen{
	public:
		typedef std::pair<int, int>                   NodeRefT;
		typedef std::vector<NodeRefT>                 RandWalkT;
		typedef boost::numeric::ublas::matrix<double> AdjMatT;

		virtual RandWalkT operator()(const AdjMatT&)=0;
};


#endif /* __RAND_WALK_GEN_HPP__ */
