
#ifndef INCLUDE__CARTESIAN_HXX

#define INCLUDE__CARTESIAN_HXX

#include <cmath>

#include <functional>
#include <numeric>

typedef double coord;
typedef std::vector<coord> point;

coord sqdiff(coord a, coord b)
{
	return (a - b) * (a - b);
}

coord calc_dist(const point &pt0, const point &pt)
{
	return std::sqrt(std::inner_product(pt0.begin(), pt0.end(), pt.begin(), 0.0,
			std::plus<coord>(), sqdiff));
}

#endif

