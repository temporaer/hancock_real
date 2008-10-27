
#include <boost/shared_ptr.hpp>
#include <boost/numeric/ublas/fwd.hpp>


class RandomAdjMatGen{
	private:
		int   mMatrixSize;
		float mConnectionProb;
		bool  mWeightedEdges;
		float mSeed;
	public:
		void configure();
		inline void setMatrixSize(int n){mMatrixSize=n;}
		inline void setConnectionProb(float f){mConnectionProb=f;}
		inline void setWeightedEdges(bool b){mWeightedEdges=b;}
		inline void setSeed(float s){mSeed=s;}

		typedef boost::numeric::ublas::matrix<double,
				boost::numeric::ublas::column_major
			> AdjMatT;
		boost::shared_ptr<AdjMatT> operator()();
};

