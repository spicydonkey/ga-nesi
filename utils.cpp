#include "utils.h"
#include <clocale>
#include <locale>
#include <vector>
#include <stdlib.h>

using namespace std;


std::string convert(const std::wstring& wstr)
{
    std::locale const loc("");
    wchar_t const *from=wstr.c_str();
    std::size_t len=wstr.size();
    std::vector<char> buffer(len+1);

    std::use_facet<std::ctype<wchar_t> >(loc).narrow(from,from+len,'_',&buffer[0]);
    return std::string(&buffer[0],&buffer[len]);    
}

std::wstring convert(const std::string& str)
{
    //wstring wstr(str.begin(),str.end());
    std::wstring wstr(str.size(), L' '); 
    wstr.resize(mbstowcs(&wstr[0], str.c_str(), str.size()));

    return wstr;
}


double rnd_generate(double min, double max)
{
    double r = (double)rand() / (double)RAND_MAX;
    return min + r * (max - min);
}

