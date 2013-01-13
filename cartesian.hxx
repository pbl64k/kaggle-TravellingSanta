
#ifndef INCLUDE__CARTESIAN_HXX

#define INCLUDE__CARTESIAN_HXX

#include <cmath>

#include <functional>
#include <numeric>
#include <utility>

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

bool xsect2(coord x1, coord y1, coord x2, coord y2, coord x3, coord y3,
		coord x4, coord y4)
{
	coord d = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);

	if (d == 0)
	{
		return false;
	}

	coord px = ((x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4)) / d;
	coord py = ((x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4)) / d;

	if ((px >= std::min(x1, x2)) && (px <= std::max(x1, x2)) &&
			(px >= std::min(x3, x4)) && (px <= std::max(x3, x4)) &&
			(py >= std::min(y1, y2)) && (py <= std::max(y1, y2)) &&
			(py >= std::min(y3, y4)) && (py <= std::max(y3, y4)))
	{
		return true;
	}

	return false;
}

#endif

