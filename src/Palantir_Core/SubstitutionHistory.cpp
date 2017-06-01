#include "SubstitutionHistory.hpp"

Palantir::SubstitutionHistory::SubstitutionHistory(): size(0) {}

Palantir::SubstitutionHistory::SubstitutionHistory(
        const deque<double> &time,
        const deque<ullong> &states_from,
        const deque<ullong> &states_to) :
        time(time),
        state_from(states_from),
        state_to(states_to),
        size(time.size()) {}

void Palantir::SubstitutionHistory::append(SubstitutionHistory &s) {
    size += s.size;
    time.insert(time.end(), s.time.begin(), s.time.end());
    state_from.insert(state_from.end(), s.state_from.begin(), s.state_from.end());
    state_to.insert(state_to.end(), s.state_to.begin(), s.state_to.end());
}


void Palantir::SubstitutionHistory::fast_forward(double t) {
    for(ullong i = 0; i < size; i++) {
        time[i] += t;
    }
}
