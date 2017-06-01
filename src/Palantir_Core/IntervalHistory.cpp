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
Palantir::IntervalHistory::IntervalHistory(
        const double end_time,
        const double length) {

    double div = floor(end_time / length);
    size = ((end_time - (length * div)) == 0) ? (ullong)div : (ullong)(div + 1);

    time_from.resize(size);
    time_to.resize(size);
    state.resize(size);

    double t = 0;
    ullong i = 0;
    while(t + length < end_time) {
        time_from[i] = t;
        time_to[i] = t + length;
        state[i] = i;

        t+=length;
        i++;
    }
    time_from.back() = t;
    time_to.back() = end_time;
    state[i] = i;
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
