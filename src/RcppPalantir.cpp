#include "RcppPalantir.hpp"

std::string get_genetic_code_name() {
    Environment globals = Environment::namespace_env("PalantiR").get(".globals");
    return globals["genetic_code_name"];
}

std::string get_attr(List l, std::string a) {
    return l[a];
}

bool has_class(List a, std::string cl) {
    std::string c = a.attr("class");
    return cl == c;
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
// set{substitution_models} =? set{intervals.state}
bool compare_modes(const List& substitution_models, const List& intervals)
{
    set<ullong> model_modes;
    for(ullong i = 0; i < substitution_models.size(); i++) {
        model_modes.insert(i);
    }

    vector<ullong> tree_states = intervals["state"];
    set<ullong> tree_modes(tree_states.begin(), tree_states.end());

    return model_modes == tree_modes;
}
