#ifndef IntervalHistory_hpp
#define IntervalHistory_hpp

#include "Util.hpp"
#include "SubstitutionHistory.hpp"

namespace Palantir
{
    class IntervalHistory
    {
        public:
            vector<ullong> state;
            vector<double> time_from;
            vector<double> time_to;
            ullong size;

            IntervalHistory();

            IntervalHistory(
                    const vector<ullong>& state,
                    const vector<double>& time_from,
                    const vector<double>& time_to);

            IntervalHistory(
                    const ullong state,
                    const double start_time,
                    const double end_time);

            IntervalHistory(
                    const SubstitutionHistory &s,
                    const ullong start_state,
                    const double end_time);

            IntervalHistory(
                    const double length,
                    const double end_time);

            void fast_forward(double time);

            void append(IntervalHistory &i);
            void append(ullong s, double from, double to);

            friend ostream& operator<<(ostream& os, IntervalHistory &it) {
                os << "state\tfrom\tto" << endl;
                for (ullong i = 0; i < it.size; i++) {
                    os << it.state[i] << "\t" << it.time_from[i] << "\t" << it.time_to[i] << endl;
                }
                return os;
            }
    };
}

#endif //IntervalHistory_hpp
