#include <cstring>
#include <iostream>
#include <sstream>
#include "utils.hpp"



void Utils::validateArg(bool condition, const std::string &msg) {
    if (!condition) {
        std::cerr << msg << std::endl;
        std::exit(1);
    }
}


std::string Utils::repeatJoin(
    const std::string &u, const std::string &s, int n)
{
    std::ostringstream out_str;
    int last = n - 1;
    for (int i = 0; i < last; ++i) {
        out_str << u << s;
    }
    out_str << u;
    return out_str.str();
}


void Utils::krepeatPuts(kstring_t *o, const std::string &u,
    const std::string &s, int n)
{
    int last = n - 1;
    for (int i = 0; i < last; ++i) {
        kputs(u.c_str(), o);
        kputs(s.c_str(), o);
    }
    kputs(u.c_str(), o);
}
