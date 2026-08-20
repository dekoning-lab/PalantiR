#include "IntervalHistory.hpp"

Palantir::IntervalHistory::IntervalHistory(): size(0) {};

Palantir::IntervalHistory::IntervalHistory(
       const deque<ullong>& state,
       const deque<double>& time_from,
       const deque<double>& time_to):
    state(state),
    time_from(time_from),
    time_to(time_to),
    size(state.size()) {}

Palantir::IntervalHistory::IntervalHistory(
        const ullong state,
        const double start_time,
        const double end_time):
    state(1, state),
    time_from(1, start_time),
    time_to(1, end_time),
    size(1) {}

Palantir::IntervalHistory::IntervalHistory(
        const SubstitutionHistory &s,
        const ullong start_state,
        const double end_time)
        :state(s.size + 1),
         time_from(s.size + 1),
         time_to(s.size + 1),
         size(s.size + 1) {

    if (s.size == 0) {
        state[0] = start_state;
        time_from[0] = 0;
        time_to[0] = end_time;
    } else {
        state[0] = s.state_from[0];
        time_from[0] = 0;
        time_to[0] = s.time[0];

        for(ullong i = 1; i < s.size; i++) {
            state[i] = s.state_from[i];
            time_from[i] = s.time[i-1];
            time_to[i] = s.time[i];
        }

        ullong last = s.size - 1;
        if (end_time != s.time[last]) {
            state[s.size] = s.state_to[last];
            time_from[s.size] = s.time[last];
            time_to[s.size] = end_time;
        }
    }
}

// Here states are ignored unless set explicitly
// FIX (2026-08-20): the closed-form `size` and the accumulating fill loop
// disagreed with each other. On (end_time, length) pairs that are exact
// multiples in double arithmetic the accumulated `t` drifts just below
// (size - 1) * length, so the loop ran `size` times instead of `size - 1`:
// `state[size]` was written one past the end of the deque, `time_from.back()`
// was overwritten with the drifted end point, and the final segment collapsed
// onto end_time -- one whole segment of branch time simply disappeared
// (end_time 0.8 at length 0.1 covered nothing in [0.7, 0.8), i.e. 12.5% of the
// branch was never simulated). The fill is now pure index arithmetic against
// the same closed-form size, so the segments tile [0, end_time] exactly with
// no accumulation drift and no out-of-range write.
// FIX (2026-08-20): degenerate inputs (end_time <= 0, which is what a
// zero-length branch hands the segment rescaler, or a non-positive length)
// used to build a size-0 history and then write to its empty deques, which
// segfaulted. They now give an empty history that the caller can skip.
Palantir::IntervalHistory::IntervalHistory(
        const double end_time,
        const double length): size(0) {

    if (!(end_time > 0) || !(length > 0)) {
        return;
    }

    double div = std::floor(end_time / length);
    size = ((end_time - (length * div)) == 0) ? (ullong)div : (ullong)(div + 1);

    time_from.resize(size);
    time_to.resize(size);
    state.resize(size);

    for (ullong i = 0; i < size; i++) {
        time_from[i] = std::min((double)i * length, end_time);
        time_to[i] = std::min((double)(i + 1) * length, end_time);
        state[i] = i;
    }
    // the last segment ends exactly at end_time whatever the rounding did
    time_to[size - 1] = end_time;
}

void Palantir::IntervalHistory::fast_forward(double time) {
    for(ullong i = 0; i < size; i++) {
        time_from[i] += time;
        time_to[i] += time;
    }
}

void Palantir::IntervalHistory::append(IntervalHistory &i) {
    size += i.size;
    state.insert(state.end(), i.state.begin(), i.state.end());
    time_from.insert(time_from.end(), i.time_from.begin(), i.time_from.end());
    time_to.insert(time_to.end(), i.time_to.begin(), i.time_to.end());
}

void Palantir::IntervalHistory::append(ullong s, double from, double to) {
    size ++;
    state.push_back(s);
    time_from.push_back(from);
    time_to.push_back(to);
}
