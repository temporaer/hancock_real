/*       Created   :  10/05/2008 07:45:51 PM
 *       Last Change: Wed Nov 12 07:00 PM 2008 CET
 */

#ifndef __SERIATION_GEN_HPP__
#define __SERIATION_GEN_HPP__
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/shared_ptr.hpp>

class SeriationGen{
	public:
		typedef std::pair<int, int>                                NodeRefT;
		typedef std::vector<NodeRefT>                              SeriationT;
		typedef boost::numeric::ublas::matrix<double,
				boost::numeric::ublas::column_major>               AdjMatT;

		virtual SeriationT operator()(boost::shared_ptr<AdjMatT>);
		virtual void configure();
};


#endif /* __SERIATION_GEN_HPP__ */
