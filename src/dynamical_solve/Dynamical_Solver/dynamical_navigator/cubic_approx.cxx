#include "../../../sdp_solve.hxx"
#include "../../../dynamical_solve.hxx"
#include <El.hpp>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <limits>
#include <vector>
//#include <gmpxx.h>
//#include <boost/multiprecision/gmp.hpp>
//#include <mpfr.h>
//#include <mpf2mpfr.h>
//using namespace boost::multiprecision;


El::BigFloat El_Pow(const El::BigFloat & a, int p)
{
	if (p > 1) return a * El_Pow(a, p - 1);
	if (p == 1) return a;
	if (p == 0) return 1;
	throw std::runtime_error(std::string("El_Pow failed: negative p"));
}

El::BigFloat CubicApprox_Optimal_x_SqrtRadicand(const El::BigFloat & x1, const El::BigFloat & y1,
	const El::BigFloat & x2, const El::BigFloat & y2,
	const El::BigFloat & x3, const El::BigFloat & y3,
	const El::BigFloat & x4, const El::BigFloat & y4)
{
	return El_Pow(x3*El_Pow(x4, 3)*(-y1 + y2) + x1 * (x1 - x4)*x4*(x1 + x4)*(y2 - y3) + El_Pow(x3, 3)*(x4*(y1 - y2) + x1 * (y2 - y4)) + El_Pow(x1, 3)*x3*(-y2 + y4) + El_Pow(x2, 3)*(-(x4*y1) - x1 * y3 + x4 * y3 + x3 * (y1 - y4) + x1 * y4) + x2 * (El_Pow(x4, 3)*(y1 - y3) + El_Pow(x1, 3)*(y3 - y4) + El_Pow(x3, 3)*(-y1 + y4)), 2) - 3 * (x3*El_Pow(x4, 2)*(-y1 + y2) + x1 * (x1 - x4)*x4*(y2 - y3) + El_Pow(x3, 2)*(x4*(y1 - y2) + x1 * (y2 - y4)) + El_Pow(x1, 2)*x3*(-y2 + y4) + El_Pow(x2, 2)*(-(x4*y1) - x1 * y3 + x4 * y3 + x3 * (y1 - y4) + x1 * y4) + x2 * (El_Pow(x4, 2)*(y1 - y3) + El_Pow(x1, 2)*(y3 - y4) + El_Pow(x3, 2)*(-y1 + y4)))*(El_Pow(x1, 2)*(x1 - x4)*El_Pow(x4, 2)*(y2 - y3) + El_Pow(x3, 3)*(El_Pow(x4, 2)*(y1 - y2) + El_Pow(x1, 2)*(y2 - y4)) + El_Pow(x2, 2)*(El_Pow(x4, 3)*(y1 - y3) + El_Pow(x1, 3)*(y3 - y4) + El_Pow(x3, 3)*(-y1 + y4)) + El_Pow(x3, 2)*(El_Pow(x4, 3)*(-y1 + y2) + El_Pow(x1, 3)*(-y2 + y4)) + El_Pow(x2, 3)*(El_Pow(x4, 2)*(-y1 + y3) + El_Pow(x3, 2)*(y1 - y4) + El_Pow(x1, 2)*(-y3 + y4)));
}


El::BigFloat CubicApprox_Optimal_x1(const El::BigFloat & x1, const El::BigFloat & y1, 
	const El::BigFloat & x2, const El::BigFloat & y2, 
	const El::BigFloat & x3, const El::BigFloat & y3, 
	const El::BigFloat & x4, const El::BigFloat & y4)
{
	return (-Sqrt(CubicApprox_Optimal_x_SqrtRadicand(x1, y1, x2, y2, x3, y3, x4, y4)) + El_Pow(x3, 3)*x4*y1 - x3 * El_Pow(x4, 3)*y1 - El_Pow(x1, 3)*x3*y2 + x1 * El_Pow(x3, 3)*y2 + El_Pow(x1, 3)*x4*y2 - El_Pow(x3, 3)*x4*y2 - x1 * El_Pow(x4, 3)*y2 + x3 * El_Pow(x4, 3)*y2 - El_Pow(x1, 3)*x4*y3 + x1 * El_Pow(x4, 3)*y3 + El_Pow(x1, 3)*x3*y4 - x1 * El_Pow(x3, 3)*y4 + x2 * (El_Pow(x4, 3)*(y1 - y3) + El_Pow(x1, 3)*(y3 - y4) + El_Pow(x3, 3)*(-y1 + y4)) + El_Pow(x2, 3)*(x4*(-y1 + y3) + x3 * (y1 - y4) + x1 * (-y3 + y4))) / (3*(x3*El_Pow(x4, 2)*(-y1 + y2) + x1 * (x1 - x4)*x4*(y2 - y3) + El_Pow(x3, 2)*(x4*(y1 - y2) + x1 * (y2 - y4)) + El_Pow(x1, 2)*x3*(-y2 + y4) + El_Pow(x2, 2)*(-(x4*y1) - x1 * y3 + x4 * y3 + x3 * (y1 - y4) + x1 * y4) + x2 * (El_Pow(x4, 2)*(y1 - y3) + El_Pow(x1, 2)*(y3 - y4) + El_Pow(x3, 2)*(-y1 + y4))));
}

El::BigFloat CubicApprox_Optimal_x2(const El::BigFloat & x1, const El::BigFloat & y1,
	const El::BigFloat & x2, const El::BigFloat & y2,
	const El::BigFloat & x3, const El::BigFloat & y3,
	const El::BigFloat & x4, const El::BigFloat & y4)
{
	return (Sqrt(CubicApprox_Optimal_x_SqrtRadicand(x1, y1, x2, y2, x3, y3, x4, y4)) + El_Pow(x3, 3)*x4*y1 - x3 * El_Pow(x4, 3)*y1 - El_Pow(x1, 3)*x3*y2 + x1 * El_Pow(x3, 3)*y2 + El_Pow(x1, 3)*x4*y2 - El_Pow(x3, 3)*x4*y2 - x1 * El_Pow(x4, 3)*y2 + x3 * El_Pow(x4, 3)*y2 - El_Pow(x1, 3)*x4*y3 + x1 * El_Pow(x4, 3)*y3 + El_Pow(x1, 3)*x3*y4 - x1 * El_Pow(x3, 3)*y4 + x2 * (El_Pow(x4, 3)*(y1 - y3) + El_Pow(x1, 3)*(y3 - y4) + El_Pow(x3, 3)*(-y1 + y4)) + El_Pow(x2, 3)*(x4*(-y1 + y3) + x3 * (y1 - y4) + x1 * (-y3 + y4))) / (3*(x3*El_Pow(x4, 2)*(-y1 + y2) + x1 * (x1 - x4)*x4*(y2 - y3) + El_Pow(x3, 2)*(x4*(y1 - y2) + x1 * (y2 - y4)) + El_Pow(x1, 2)*x3*(-y2 + y4) + El_Pow(x2, 2)*(-(x4*y1) - x1 * y3 + x4 * y3 + x3 * (y1 - y4) + x1 * y4) + x2 * (El_Pow(x4, 2)*(y1 - y3) + El_Pow(x1, 2)*(y3 - y4) + El_Pow(x3, 2)*(-y1 + y4))));
}

El::BigFloat CubicApprox_Optimal_x(const El::BigFloat & SqrtPart, const El::BigFloat & x1, const El::BigFloat & y1,
	const El::BigFloat & x2, const El::BigFloat & y2,
	const El::BigFloat & x3, const El::BigFloat & y3,
	const El::BigFloat & x4, const El::BigFloat & y4)
{
	return (SqrtPart + El_Pow(x3, 3)*x4*y1 - x3 * El_Pow(x4, 3)*y1 - El_Pow(x1, 3)*x3*y2 + x1 * El_Pow(x3, 3)*y2 + El_Pow(x1, 3)*x4*y2 - El_Pow(x3, 3)*x4*y2 - x1 * El_Pow(x4, 3)*y2 + x3 * El_Pow(x4, 3)*y2 - El_Pow(x1, 3)*x4*y3 + x1 * El_Pow(x4, 3)*y3 + El_Pow(x1, 3)*x3*y4 - x1 * El_Pow(x3, 3)*y4 + x2 * (El_Pow(x4, 3)*(y1 - y3) + El_Pow(x1, 3)*(y3 - y4) + El_Pow(x3, 3)*(-y1 + y4)) + El_Pow(x2, 3)*(x4*(-y1 + y3) + x3 * (y1 - y4) + x1 * (-y3 + y4))) / (3*(x3*El_Pow(x4, 2)*(-y1 + y2) + x1 * (x1 - x4)*x4*(y2 - y3) + El_Pow(x3, 2)*(x4*(y1 - y2) + x1 * (y2 - y4)) + El_Pow(x1, 2)*x3*(-y2 + y4) + El_Pow(x2, 2)*(-(x4*y1) - x1 * y3 + x4 * y3 + x3 * (y1 - y4) + x1 * y4) + x2 * (El_Pow(x4, 2)*(y1 - y3) + El_Pow(x1, 2)*(y3 - y4) + El_Pow(x3, 2)*(-y1 + y4))));
}


El::BigFloat CubicApprox_Func(const El::BigFloat & x1, const El::BigFloat & y1,
	const El::BigFloat & x2, const El::BigFloat & y2,
	const El::BigFloat & x3, const El::BigFloat & y3,
	const El::BigFloat & x4, const El::BigFloat & y4, const El::BigFloat & x)
{
	return (x1*(x1 - x3)*x3*(x1 - x4)*(x3 - x4)*x4*y2 + El_Pow(x2, 3)*(x1*(x1 - x4)*x4*y3 + El_Pow(x3, 2)*(-(x4*y1) + x1 * y4) + x3 * (El_Pow(x4, 2)*y1 - El_Pow(x1, 2)*y4)) + x2 * (El_Pow(x1, 2)*(x1 - x4)*El_Pow(x4, 2)*y3 + El_Pow(x3, 3)*(-(El_Pow(x4, 2)*y1) + El_Pow(x1, 2)*y4) + El_Pow(x3, 2)*(El_Pow(x4, 3)*y1 - El_Pow(x1, 3)*y4)) + El_Pow(x2, 2)*(x1*x4*(-El_Pow(x1, 2) + El_Pow(x4, 2))*y3 + El_Pow(x3, 3)*(x4*y1 - x1 * y4) + x3 * (-(El_Pow(x4, 3)*y1) + El_Pow(x1, 3)*y4)) + El_Pow(x, 3)*(x1*(x1 - x4)*x4*(y2 - y3) + El_Pow(x3, 2)*(x4*(y1 - y2) + x1 * (y2 - y4)) + x2 * (El_Pow(x4, 2)*(y1 - y3) + El_Pow(x1, 2)*(y3 - y4) + El_Pow(x3, 2)*(-y1 + y4)) + x3 * (El_Pow(x4, 2)*(-y1 + y2) + El_Pow(x1, 2)*(-y2 + y4)) + El_Pow(x2, 2)*(x4*(-y1 + y3) + x3 * (y1 - y4) + x1 * (-y3 + y4))) + x * (El_Pow(x1, 2)*(x1 - x4)*El_Pow(x4, 2)*(y2 - y3) + El_Pow(x3, 3)*(El_Pow(x4, 2)*(y1 - y2) + El_Pow(x1, 2)*(y2 - y4)) + El_Pow(x2, 2)*(El_Pow(x4, 3)*(y1 - y3) + El_Pow(x1, 3)*(y3 - y4) + El_Pow(x3, 3)*(-y1 + y4)) + El_Pow(x3, 2)*(El_Pow(x4, 3)*(-y1 + y2) + El_Pow(x1, 3)*(-y2 + y4)) + El_Pow(x2, 3)*(El_Pow(x4, 2)*(-y1 + y3) + El_Pow(x3, 2)*(y1 - y4) + El_Pow(x1, 2)*(-y3 + y4))) + El_Pow(x, 2)*(-(x1*x4*(El_Pow(x1, 2) - El_Pow(x4, 2))*(y2 - y3)) + x3 * (El_Pow(x4, 3)*(y1 - y2) + El_Pow(x1, 3)*(y2 - y4)) + El_Pow(x2, 3)*(x4*(y1 - y3) + x1 * (y3 - y4) + x3 * (-y1 + y4)) + El_Pow(x3, 3)*(x4*(-y1 + y2) + x1 * (-y2 + y4)) + x2 * (El_Pow(x4, 3)*(-y1 + y3) + El_Pow(x3, 3)*(y1 - y4) + El_Pow(x1, 3)*(-y3 + y4)))) / ((x1 - x2)*(x1 - x3)*(x2 - x3)*(x1 - x4)*(x2 - x4)*(x3 - x4));
}

bool CubicApprox_Max_x(const El::BigFloat & x1, const El::BigFloat & y1,
	const El::BigFloat & x2, const El::BigFloat & y2,
	const El::BigFloat & x3, const El::BigFloat & y3,
	const El::BigFloat & x4, const El::BigFloat & y4,
	El::BigFloat & max_x)
{
	El::BigFloat optimal_radicand = CubicApprox_Optimal_x_SqrtRadicand(x1, y1, x2, y2, x3, y3, x4, y4);
	if (optimal_radicand < 0)
	{
		std::vector<El::BigFloat> ylist{ y1,y2,y3,y4 };
		std::vector<El::BigFloat> xlist{ x1,x2,x3,x4 };

		auto max_y_it = max_element(ylist.begin(), ylist.end());
		int max_index = max_y_it - ylist.begin();
		max_x = xlist.at(max_index);
		return false;
	}

	El::BigFloat optimal_sqrtpart = El::Sqrt(optimal_radicand);

	El::BigFloat optimal_x1 = CubicApprox_Optimal_x(optimal_sqrtpart, x1, y1, x2, y2, x3, y3, x4, y4);
	El::BigFloat optimal_x2 = CubicApprox_Optimal_x(-optimal_sqrtpart, x1, y1, x2, y2, x3, y3, x4, y4);

	El::BigFloat optimal_y1 = CubicApprox_Func(x1, y1, x2, y2, x3, y3, x4, y4, optimal_x1);
	El::BigFloat optimal_y2 = CubicApprox_Func(x1, y1, x2, y2, x3, y3, x4, y4, optimal_x2);

	max_x = (optimal_y1 > optimal_y2) ? optimal_x1 : optimal_x2;

	return true;
}


template <class T1, class T2, class Pred = std::less<T2> >
struct comp_pair_second {
	bool operator()(const std::pair<T1, T2>&left, const std::pair<T1, T2>&right) {
		Pred p;
		return p(left.second, right.second);
	}
};
template <class T1, class T2, class Pred = std::less<T1> >
struct comp_pair_first {
	bool operator()(const std::pair<T1, T2>&left, const std::pair<T1, T2>&right) {
		Pred p;
		return p(left.first, right.first);
	}
};

bool CubicApprox_LineSearch_RecommendPoint(std::vector<std::pair<El::BigFloat, El::BigFloat>> & history, 
	const El::BigFloat & beta_scan_end, const El::BigFloat & beta_scan_step, El::BigFloat & new_x)
{
	if (history.size() > 20)return false;

	El::BigFloat x_max = max_element(history.begin(), history.end(), comp_pair_first<El::BigFloat, El::BigFloat>())->first;

	if (x_max <= beta_scan_end && x_max + beta_scan_step > beta_scan_end)
	{
		new_x = x_max + beta_scan_step;
		return true;
	}

	if (max_element(history.begin(), history.end(), comp_pair_second<El::BigFloat, El::BigFloat>())->second > 0.7) return false;

	bool cubic_approx_validQ = CubicApprox_Max_x(history.end()[-4].first, history.end()[-4].second,
		history.end()[-3].first, history.end()[-3].second,
		history.end()[-2].first, history.end()[-2].second,
		history.end()[-1].first, history.end()[-1].second,
		new_x);

	if(cubic_approx_validQ ==false) return false;

	El::BigFloat range_min = min_element(history.end() - 4, history.end(), comp_pair_first<El::BigFloat, El::BigFloat>())->first;
	El::BigFloat range_max = max_element(history.end() - 4, history.end(), comp_pair_first<El::BigFloat, El::BigFloat>())->first;

	if (new_x<range_max && new_x>range_min)return true;

	return false;
}

