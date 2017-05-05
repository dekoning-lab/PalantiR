#include "RcppPalantir.hpp"

std::string get_genetic_code_name() {
    Environment globals = Environment::namespace_env("PalantiR").get(".globals");
    return globals["genetic_code_name"];
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
