
#ifndef INCLUDE__QXT_HXX

#define INCLUDE__QXT_HXX

#include <cassert>

#include <map>
#include <set>
#include <utility>

#include "cartesian.hxx"

typedef std::pair<point, point> qx_edge;

template<typename T> class qxt
{
	public:
	bool top_;
	size_t depth_;
	point tl_, br_, mid_;
	std::map<qx_edge, T> es_;
	qxt<T> *q1_;
	qxt<T> *q2_;
	qxt<T> *q3_;
	qxt<T> *q4_;

	qxt():
			top_(true), depth_(10), tl_(2), br_(2), mid_(),
			es_(), q1_(NULL), q2_(NULL), q3_(NULL), q4_(NULL)
	{
		tl_[0] = 0.0;
		tl_[1] = 0.0;
		br_[0] = 20000.0;
		br_[1] = 20000.0;
		
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

	qxt(bool top, size_t depth, point tl, point br):
			top_(top), depth_(depth), tl_(tl), br_(br), mid_(),
			es_(), q1_(NULL), q2_(NULL), q3_(NULL), q4_(NULL)
	{
#ifndef NDEBUG
		assert(depth_ < 24);
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

	qxt(const qxt<T> &orig): top_(orig.top_), depth_(orig.depth_),
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
			top_ = rhs.top_;
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

	void add_edge(const qx_edge &e, const T& data)
	{
#ifndef NDEBUG
		assert(top_);
#endif

		if ((e.first[0] > e.second[0]) || ((e.first[0] == e.second[0]) && (e.first[1] > e.second[1])))
		{
			add_edge0(make_pair(e.second, e.first), data);
		}
		else
		{
			add_edge0(e, data);
		}
	}

	void add_edge0(const qx_edge &e, const T& data)
	{
#ifndef NDEBUG
		if (top_)
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

		if (intersects(e))
		{
			if (depth_ == 0)
			{
				es_[e] = data;
			}
			else
			{
				q1_->add_edge0(e, data);
				q2_->add_edge0(e, data);
				q3_->add_edge0(e, data);
				q4_->add_edge0(e, data);
			}
		}
	}

	void del_edge(const qx_edge &e)
	{
#ifndef NDEBUG
		assert(top_);
#endif

		if ((e.first[0] > e.second[0]) || ((e.first[0] == e.second[0]) && (e.first[1] > e.second[1])))
		{
			del_edge0(make_pair(e.second, e.first));
		}
		else
		{
			del_edge0(e);
		}
	}

	void del_edge0(const qx_edge &e)
	{
#ifndef NDEBUG
		if (top_)
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

		if (intersects(e))
		{
			if (depth_ == 0)
			{
				assert(es_.erase(e) == 1);
			}
			else
			{
				q1_->del_edge0(e);
				q2_->del_edge0(e);
				q3_->del_edge0(e);
				q4_->del_edge0(e);
			}
		}
	}

	std::set<T> find_x(const qx_edge &e) const
	{
#ifndef NDEBUG
		assert(top_);

		assert(e.first[0] >= tl_[0]);
		assert(e.first[1] >= tl_[1]);
		assert(e.first[0] <= br_[0]);
		assert(e.first[1] <= br_[1]);
		assert(e.second[0] >= tl_[0]);
		assert(e.second[1] >= tl_[1]);
		assert(e.second[0] <= br_[0]);
		assert(e.second[1] <= br_[1]);
#endif

		std::set<T> res;

		find_x0(e, res);

		return res;
	}

	void find_x0(const qx_edge &e, std::set<T> &res) const
	{
		if (intersects(e))
		{
			if (depth_ == 0)
			{
				for (typename std::map<qx_edge, T>::const_iterator iter = es_.begin();
						iter != es_.end(); ++iter)
				{
					if (xsect2(e.first[0], e.first[1], e.second[0], e.second[1],
							iter->first.first[0], iter->first.first[1], iter->first.second[0], iter->first.second[1]))
					{
						res.insert(iter->second);
					}
				}
			}
			else
			{
				q1_->find_x0(e, res);
				q2_->find_x0(e, res);
				q3_->find_x0(e, res);
				q4_->find_x0(e, res);
			}
		}
	}

	bool intersects(const qx_edge &e) const
	{
		if (((e.first[0] >= tl_[0]) &&
				(e.first[0] <= br_[0]) &&
				(e.first[1] >= tl_[1]) &&
				(e.first[1] <= br_[1])) ||
				((e.second[0] >= tl_[0]) &&
				(e.second[0] <= br_[0]) &&
				(e.second[1] >= tl_[1]) &&
				(e.second[1] <= br_[1])))
		{
			return true;
		}

		return xsect2(e.first[0], e.first[1], e.second[0], e.second[1], tl_[0], tl_[1], br_[0], tl_[1]) ||
				xsect2(e.first[0], e.first[1], e.second[0], e.second[1], tl_[0], tl_[1], tl_[0], br_[1]) ||
				xsect2(e.first[0], e.first[1], e.second[0], e.second[1], br_[0], tl_[1], br_[0], br_[1]) ||
				xsect2(e.first[0], e.first[1], e.second[0], e.second[1], tl_[0], br_[1], br_[0], br_[1]);
	}
};

#endif

