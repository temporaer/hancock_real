/*       Created   :  10/05/2008 10:21:12 PM
 *       Last Change: Sun Oct 05 10:00 PM 2008 CEST
 */

#ifndef __SDP_RAND_WALK_GEN_HPP__
#define __SDP_RAND_WALK_GEN_HPP__
#include <rand_walk_gen.hpp>
class SDPRandWalkGen : public RandWalkGen 
{
	public:
		virtual RandWalkT operator()(const AdjMatT&);
};
#endif /* __SDP_RAND_WALK_GEN_HPP__ */

