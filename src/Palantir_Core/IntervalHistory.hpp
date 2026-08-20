#ifndef IntervalHistory_hpp
#define IntervalHistory_hpp

#include "Util.hpp"
#include "SubstitutionHistory.hpp"

namespace Palantir
{
    class IntervalHistory
    {
        public:
            deque<ullong> state;
            deque<double> time_from;
            deque<double> time_to;
            ullong size;

            IntervalHistory();

            IntervalHistory(
                    const deque<ullong>& state,
                    const deque<double>& time_from,
                    const deque<double>& time_to);

            IntervalHistory(
                    const ullong state,
                    const double start_time,
                    const double end_time);

            IntervalHistory(
                    const SubstitutionHistory &s,
                    const ullong start_state,
                    const double end_time);

            // FIX (2026-08-20): the declaration named these (length, end_time)
            // while the definition reads them as (end_time, length). Callers
            // that trusted the header passed them the wrong way round.
            IntervalHistory(
                    const double end_time,
                    const double length);

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
