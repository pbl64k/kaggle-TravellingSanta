
#ifndef INCLUDE__QXT_HXX

#define INCLUDE__QXT_HXX

#include <cassert>

#include <list>
#include <utility>

#include "cartesian.hxx"

typedef std::pair<point, point> qx_edge;

template<typename T> class qxt
{
	public:
	bool sanity_;
	size_t depth_;
	point tl_, br_, mid_;
	std::list<std::pair<qx_edge, T> > es_;
	qxt<T> *q1_;
	qxt<T> *q2_;
	qxt<T> *q3_;
	qxt<T> *q4_;

	qxt(bool sanity, size_t depth, point tl, point br):
			sanity_(sanity), depth_(depth), tl_(tl), br_(br), mid_(),
			es_(), q1_(NULL), q2_(NULL), q3_(NULL), q4_(NULL)
	{
#ifndef NDEBUG
		assert(depth_ >= 0);
		assert(tl_.size() == 2);
		assert(br_.size() == 2);
		assert(tl_[0] < br_[0]);
		assert(tl_[1] < br_[1]);
#endif
		
		if (depth_ > 0)
		{
			coord xm = (tl_[0] / 2.0) + (br_[0] / 2.0);
			coord ym = (tl_[1] / 2.0) + (br_[1] / 2.0);
			point mid(2);
			mid[0] = xm;
			mid[1] = ym;

			mid_ = mid;

			q1_ = new qxt<T>(false, depth_ - 1, tl_, mid_);
			q4_ = new qxt<T>(false, depth_ - 1, mid_, br_);
			point p1(2);
			point p2(2);
			p1[0] = mid_[0];
			p1[1] = tl_[1];
			p2[0] = br_[0];
			p2[1] = mid_[1];
			q2_ = new qxt<T>(false, depth_ - 1, p1, p2);
			p1[0] = tl_[0];
			p1[1] = mid_[1];
			p2[0] = mid_[0];
			p2[1] = br_[1];
			q3_ = new qxt<T>(false, depth_ - 1, p1, p2);
		}
	}

	qxt(const qxt<T> &orig): sanity_(orig.sanity_), depth_(orig.depth_),
			tl_(orig.tl_), br_(orig.br_), mid_(orig.mid_), es_(orig.es_),
			q1_(orig.q1_), q2_(orig.q2_), q3_(orig.q3_), q4_(orig.q4_)
	{
#ifndef NDEBUG
		assert(false);
#endif
	}

	~qxt()
	{
		if (depth_ > 0)
		{
			delete q1_;
			delete q2_;
			delete q3_;
			delete q4_;
		}
	}

	qxt<T> &operator=(const qxt<T> &rhs)
	{
#ifndef NDEBUG
		assert(false);
#endif

		if (this != (&rhs))
		{
			sanity_ = rhs.sanity_;
			depth_ = rhs.depth_;
			tl_ = rhs.tl_;
			br_ = rhs.br_;
			mid_ = rhs.mid_;
			es_ = rhs.es_;
			q1_ = rhs.q1_;
			q2_ = rhs.q2_;
			q3_ = rhs.q3_;
			q4_ = rhs.q4_;
		}

		return *this;
	}

	void add_edge(const qx_edge &e)
	{
#ifndef NDEBUG
		if (sanity_)
		{
			assert(e.first[0] >= tl_[0]);
			assert(e.first[1] >= tl_[1]);
			assert(e.first[0] <= br_[0]);
			assert(e.first[1] <= br_[1]);
			assert(e.second[0] >= tl_[0]);
			assert(e.second[1] >= tl_[1]);
			assert(e.second[0] <= br_[0]);
			assert(e.second[1] <= br_[1]);
		}
#endif
	}
};

#endif

