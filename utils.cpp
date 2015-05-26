#include "utils.h"
#include <clocale>
#include <locale>
#include <vector>
#include <stdlib.h>

using namespace std;

// convert wstring to a char string with '_' as default char
std::string convert(const std::wstring& wstr)
{
    std::locale const loc("");	// construct a locale obj set to the env's default
    wchar_t const *from=wstr.c_str();	// ptr to const wide character argument
    std::size_t len=wstr.size();	// number of characters in the wide string argument
    std::vector<char> buffer(len+1);	// init char vector buffer of sufficient length to store the input wstring

	// fill buffer with wchar_t->char transformed wstr, using '_' as default character
    std::use_facet<std::ctype<wchar_t> >(loc).narrow(from,from+len,'_',&buffer[0]);

    return std::string(&buffer[0],&buffer[len]);	// return the translated string
}

// convert string to a wstring
std::wstring convert(const std::string& str)
{
    //wstring wstr(str.begin(),str.end());
    std::wstring wstr(str.size(), L' ');	// initialise a wide string filled with wchar ' ', length equal to input string
    wstr.resize(mbstowcs(&wstr[0], str.c_str(), str.size()));	// convert input string to a w-char str pointed by wstr and resize to the number of chars translated

    return wstr;	// return the translated wstring
}

// rnd_generate
// generates a random double in [min,max]
double rnd_generate(double min, double max)
{
    double r = (double)rand() / (double)RAND_MAX;
    return min + r * (max - min);
}

