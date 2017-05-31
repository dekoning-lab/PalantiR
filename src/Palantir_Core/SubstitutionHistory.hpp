#ifndef SubstitutionHistory_hpp
#define SubstitutionHistory_hpp

#include "Util.hpp"

namespace Palantir
{
    class SubstitutionHistory 
    {
        public:

            vector<double> time;
            vector<ullong> state_from;
            vector<ullong> state_to;
            ullong size;

            SubstitutionHistory();

            SubstitutionHistory(
                    const vector<double> &time,
                    const vector<ullong> &states_from,
                    const vector<ullong> &states_to);

            void append(SubstitutionHistory &s);

            void fast_forward(double t);

            friend ostream& operator<<(ostream& os, SubstitutionHistory &s) {
                os << "from\tto\ttime" << endl;
                for (ullong i = 0; i < s.size; i++) {
                    os << s.state_from[i] << "\t" << s.state_to[i] << "\t" << s.time[i] << endl;
                }
                return os;
            }
    };
}

#endif //SubstitutionHistory_hpp
