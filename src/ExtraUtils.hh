/*
 * ExtraUtils.hh
 *
 *  created on: 06.06.2017
 *      author: M.Khaled
 */

#ifndef EXTRAUTILS_HH_
#define EXTRAUTILS_HH_

#include <cmath>
#include <array>

class ExtraUtils {
public:

	/* compare two dp-floating-point vars up ti precision */
	static
	inline
	bool doubleEquals(double a, double b, double epsilon = 0.001)
	{
	    return std::fabs(a - b) < epsilon;
	}


	static
	inline
	bool doubleEquals(size_t size, const double*  a, const double* b, double epsilon = 0.001)
	{
		for(size_t i=0; i<size; i++)
			if(!doubleEquals(a[i],b[i], epsilon))
				return false;

		return true;
	}

};
#endif
