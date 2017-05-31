#ifndef Util_hpp
#define Util_hpp

#ifdef is_R
#include <RcppArmadillo.h>
#else
#include <armadillo>
#endif

#include <string>
#include <algorithm>
#include <deque>
#include <vector>
#include <map>
#include <functional>
#include <iostream>
#include <sstream>

typedef unsigned long long ullong;

#define SMALL 1e-12
#define PUNY 1e-30
#define TINY 1e-100

using namespace std;
using namespace arma;

class Match
{
    public:
        ullong begin;
        ullong length;
        Match(): begin(string::npos), length(0) {}

        Match(ullong b, ullong l): begin(b), length(l) {}

        bool operator ==(const Match& rhs)
        {
            return ((begin == rhs.begin) && (length == rhs.length));
        }

        bool operator !=(const Match& rhs)
        {
            return !((*this) == rhs);
        }

        friend ostream& operator<<(ostream& os, const Match& m)
        {
            os << "Match(begin=" << m.begin << ", length=" << m.length << ")";
            return os;
        }
};



static Match noMatch(string::npos, 0);

void removeWhitespace(string &str);

Match find_first_in(string& s, vector<string> pattern, ullong pos = 0);
vector<string> tokenize(string& str, vector<string> delimiters, bool keepDelimiters = true, bool keepEmpty = true);

string file_to_string(const string path);

vector<string> chars(const string& s);

#endif /* Util_hpp */
