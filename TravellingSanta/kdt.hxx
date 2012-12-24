
#ifndef INCLUDE__KDT_HXX

#define INCLUDE__KDT_HXX

#include <cassert>
#include <cmath>

#include <deque>
#include <numeric>
#include <utility>
#include <vector>

typedef double coord;
typedef std::vector<coord> point;

coord sqdiff(coord a, coord b)
{
	return (a - b) * (a - b);
}

coord calc_dist(const point &pt0, const point &pt)
{
	return std::sqrt(std::inner_product(pt0.begin(), pt0.end(), pt.begin(), 0.0, std::plus<coord>(), sqdiff));
}

template<size_t N, typename T> class kdt_node
{
	public:
	point pt_;
	T data_;
	kdt_node<N, T> *lc_, *rc_;

	kdt_node(const point &pt, const T &data):
			pt_(pt), data_(data), lc_(NULL), rc_(NULL)
	{
#ifndef NDEBUG
		assert(pt_.size() == N);
#endif
	}
};

template<size_t N, typename T> class kdt
{
	public:
	kdt_node<N, T> *root_;

	kdt(): root_(NULL)
	{
	}

	void add(const point &pt, const T &data)
	{
		root_ = add0(pt, data, root_, 0);
	}

	std::deque<std::pair<coord, T> > find_knn(const point &pt, size_t k, bool (*f)(const T &)) const
	{
#ifndef NDEBUG
		assert(pt.size() == N);
#endif

		std::deque<std::pair<coord, T> > res;

		find_knn0(pt, k, root_, 0, res, f);

		return res;
	}

	kdt_node<N, T> *add0(const point &pt, const T &data, kdt_node<N, T> *node, size_t ix)
	{
		if (node == NULL)
		{
			return new kdt_node<N, T>(pt, data);
		}

		if (pt[ix] < node->pt_[ix])
		{
			node->lc_ = add0(pt, data, node->lc_, (ix + 1) % N);
		}
		else
		{
			node->rc_ = add0(pt, data, node->rc_, (ix + 1) % N);
		}

		return node;
	}

	void find_knn0(const point &pt, size_t k, kdt_node<N, T> *node, size_t ix, std::deque<std::pair<coord, T> > &res,
			bool (*f)(const T &)) const
	{
		if (node == NULL)
		{
			return;
		}

		inject_nn(pt, node->pt_, k, node->data_, res, f);

		if (pt[ix] < node->pt_[ix])
		{
			find_knn0(pt, k, node->lc_, (ix + 1) % N, res, f);

			if ((res.size() < k) || (res.back().first > abs(node->pt_[ix] - pt[ix])))
			{
				find_knn0(pt, k, node->rc_, (ix + 1) % N, res, f);
			}
		}
		else
		{
			find_knn0(pt, k, node->rc_, (ix + 1) % N, res, f);

			if ((res.size() < k) || (res.back().first > abs(node->pt_[ix] - pt[ix])))
			{
				find_knn0(pt, k, node->lc_, (ix + 1) % N, res, f);
			}
		}
	}

	void inject_nn(const point &pt0, const point &pt, size_t k, const T &data,
			std::deque<std::pair<coord, T> > &res, bool (*f)(const T &)) const
	{
		bool p = (*f)(data);

		if (! p)
		{
			return;
		}

		coord dist = calc_dist(pt0, pt);

		if (res.empty())
		{
			res.push_back(std::make_pair(dist, data));

			return;
		}

		bool injected = false;

		for (typename std::deque<std::pair<coord, T> >::iterator iter = res.begin(); iter != res.end(); ++iter)
		{
			if (dist < iter->first)
			{
				res.insert(iter, std::make_pair(dist, data));

				injected = true;

				break;
			}
		}

		if (! injected)
		{
			res.push_back(std::make_pair(dist, data));
		}

		if (res.size() > k)
		{
			res.pop_back();
		}
	}
};

#endif

