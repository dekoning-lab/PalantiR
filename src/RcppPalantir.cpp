#include "RcppPalantir.hpp"

std::string get_genetic_code_name() {
    Environment globals = Environment::namespace_env("PalantiR").get(".globals");
    return globals["genetic_code_name"];
}

// FIX (2026-08-20, MIN7): `return l[a];` threw Rcpp's opaque "Object was
// created without names" (unnamed list) or "no such index" (named list, wrong
// field) whenever a caller passed the wrong kind of object, which is precisely
// when a useful message is needed. Name the missing field instead.
std::string get_attr(List l, std::string a) {
    if(l.size() == 0 || !l.containsElementNamed(a.c_str())) {
        stop("Object is missing the required field `" + a + "`");
    }
    return as<std::string>(l[a]);
}

// FIX (2026-08-20, MIN7): `std::string c = a.attr("class");` threw on any
// object with no class attribute rather than answering the question that was
// asked. It also compared against only the first class of a multi-class object.
// Return false for a classless object and test membership of the whole vector.
bool has_class(List a, std::string cl) {
    RObject klass = a.attr("class");
    if(klass.isNULL()) {
        return false;
    }
    CharacterVector classes = as<CharacterVector>(klass);
    for(R_xlen_t i = 0; i < classes.size(); i++) {
        if(cl == as<std::string>(classes[i])) {
            return true;
        }
    }
    return false;
}

// FIX (2026-08-20, M3): counts (population sizes, sequence lengths) used to be
// declared `unsigned long long` on the Rcpp boundary, so R doubles were coerced
// before anything could look at them: 0 passed straight through (an all-NaN
// model with a normal-looking print), -1 wrapped to 1.8e19, and 2.5 was
// silently truncated. Take the argument as a double and validate it here,
// BEFORE the conversion.
unsigned long long as_count(double value, std::string argument, unsigned long long minimum)
{
    if(!R_finite(value)) {
        stop("Argument `" + argument + "` should be a finite number, not " +
             std::to_string(value));
    }
    if(value != std::floor(value)) {
        stop("Argument `" + argument + "` should be a whole number, not " +
             std::to_string(value));
    }
    if(value < (double) minimum) {
        stop("Argument `" + argument + "` should be at least " +
             std::to_string(minimum) + ", not " + std::to_string((long long) value));
    }
    if(value > 9007199254740992.0) {
        stop("Argument `" + argument + "` is too large to represent exactly");
    }
    return (unsigned long long) value;
}

// FIX (2026-08-20, M3): mutation rates of 0 or below give a degenerate or
// sign-flipped rate matrix that then propagates NaN through the scaling.
void as_positive_rate(double value, std::string argument)
{
    if(!R_finite(value) || value <= 0) {
        stop("Argument `" + argument + "` should be a finite number greater "
             "than 0, not " + std::to_string(value));
    }
}

// FIX (2026-08-20, C8): armadillo's unchecked accessors (.at()/operator[]) read
// out of bounds when a fitness vector or delta matrix is the wrong size, which
// produced models with non-finite rates that differed between calls. Check the
// shape up front and name the argument, the expected size and what arrived.
void check_length(unsigned long long actual, unsigned long long expected, std::string argument)
{
    if(actual != expected) {
        stop("Argument `" + argument + "` should have " +
             std::to_string(expected) + " elements (one per amino acid under "
             "the active genetic code), not " + std::to_string(actual));
    }
}

void check_square(unsigned long long rows, unsigned long long cols,
                  unsigned long long expected, std::string argument)
{
    if(rows != expected || cols != expected) {
        stop("Argument `" + argument + "` should be a " +
             std::to_string(expected) + "x" + std::to_string(expected) +
             " matrix (one row and column per amino acid under the active "
             "genetic code), not " + std::to_string(rows) + "x" +
             std::to_string(cols));
    }
}

CharacterVector predicate(arma::uvec predicate, string yes, string no)
{
    CharacterVector c(predicate.size());
    for(unsigned long long i = 0; i < predicate.size(); i++) {
        if(predicate[i]) {
            c[i] = yes;
        } else {
            c[i] = no;
        }
    }
    return c;
}

List interval_histories_to_list(const std::vector<Palantir::IntervalHistory>& tree_intervals, const Palantir::Phylogeny& p)
{
    vector<ullong> node_index;
    vector<ullong> state;
    vector<double> time_from;
    vector<double> time_to;

    for(ullong node = 0; node < p.n_nodes; node++) {
        const Palantir::IntervalHistory& intervals = tree_intervals[node];

        state.insert(state.end(), intervals.state.begin(), intervals.state.end());
        time_from.insert(time_from.end(), intervals.time_from.begin(), intervals.time_from.end());
        time_to.insert(time_to.end(), intervals.time_to.begin(), intervals.time_to.end());

        vector<ullong> t(intervals.size, node);
        node_index.insert(node_index.end(), t.begin(), t.end());
    }

    return List::create(
        _["node"] = node_index,
        _["state"] = state,
        _["from"] = time_from,
        _["to"] = time_to);
}

List site_simulations_to_list(const std::vector<Palantir::SiteSimulation>& sims, const Palantir::Phylogeny& p)
{
    vector<ullong> node_index;
    vector<ullong> sites;
    vector<double> time;
    vector<ullong> state_from;
    vector<ullong> state_to;

    for(ullong site = 0; site < sims.size(); site++) {
        ullong n_site_subs = 0;
        for(ullong node = 0; node < p.n_nodes; node++) {
            const Palantir::SubstitutionHistory& s = sims[site].substitutions[node];

            time.insert(time.end(), s.time.begin(), s.time.end());
            state_from.insert(state_from.end(), s.state_from.begin(), s.state_from.end());
            state_to.insert(state_to.end(), s.state_to.begin(), s.state_to.end());

            vector<ullong> t(s.size, node);
            node_index.insert(node_index.end(), t.begin(), t.end());

            n_site_subs += s.size;
        }
        vector<ullong> t(n_site_subs, site);
        sites.insert(sites.end(), t.begin(), t.end());
    }

    return List::create(
        _["site"] = sites,
        _["node"] = node_index,
        _["time"] = time,
        _["from"] = state_from,
        _["to"] = state_to);
}

// Do substitution models and intervals have the same modes?
// set{substitution_models} contains set{intervals.state}
//[[Rcpp::export]]
bool compare_modes(const List& substitution_models, const List& intervals)
{
    vector<ullong> model_modes;
    for(ullong i = 0; i < substitution_models.size(); i++)
        model_modes.push_back(i);

    vector<ullong> tree_states = intervals["state"];
    set<ullong> unique_tree_states(tree_states.begin(), tree_states.end());
    return includes(model_modes.begin(), model_modes.end(),
                    unique_tree_states.begin(), unique_tree_states.end());
}
