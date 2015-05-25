#ifndef UTILS_H
#define UTILS_H
#include <string>
#include <limits>

#define MAX_DOUBLE std::numeric_limits<double>::max()


std::string convert(const std::wstring& wstr);
std::wstring convert(const std::string& str);

double rnd_generate(double min, double max);

template<class T,class S> struct pair_equal_to:std::binary_function<T,std::pair<T,S>,bool> {
  bool operator()(const T& y, const std::pair<T,S>& x) const
  {
      return x.first==y;
  }
};


//returns true if value is within +/- eps from point
inline bool in_range(double point, double value, double eps)
{
    return (point>=value-eps && point<=value+eps);
}

#endif

