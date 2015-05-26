#ifndef UTILS_H
#define UTILS_H
#include <string>
#include <limits>

#define MAX_DOUBLE std::numeric_limits<double>::max()


std::string convert(const std::wstring& wstr);	// convert wstring to a char string with '_' as default char
std::wstring convert(const std::string& str);	// convert string to a wstring

// generate a random double in [min,max]
double rnd_generate(double min, double max);


//pair_equal_to
//contains the binary operator to evaluate if a pair is equal to an obj (type of first memb)
template<class T,class S> struct pair_equal_to:std::binary_function<T,std::pair<T,S>,bool> {
	bool operator()(const T& y, const std::pair<T,S>& x) const
	{
		return x.first==y;	//pair equal to &T iff its first member is equal to &T
	}
};


//returns true if value is within +/- eps from point
inline bool in_range(double point, double value, double eps)
{
    return (point>=value-eps && point<=value+eps);
}

#endif

