#include <RcppArmadillo.h>

#include "GeneticCode.hpp"
#include "Palantir.hpp"

#include "RcppUtil.hpp"

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
List Sequence(std::string sequence, std::string type = "codon", unsigned long long mode = 0)
{
    Palantir::GeneticCode g(get_genetic_code_name());

    vector<string> s;
    uvec i;
    if(type == "nucleotide") {
        s = Palantir::Nucleotide::split(sequence);
        i = Palantir::Nucleotide::to_digit(s);
    }
    else if(type == "amino_acid") {
        s = Palantir::AminoAcid::split(sequence);
        i = Palantir::AminoAcid::to_digit(s);
    }
    else if(type == "codon") {
        s = Palantir::Codon::split(sequence);
        i = Palantir::Codon::to_digit(s, g);
    }
    else if(type == "codon_pair") {
        s = Palantir::CodonPair::split(sequence);
        i = Palantir::CodonPair::to_digit(s, g);
    }
    else if(type == "compound_codon") {
        s = Palantir::CompoundCodon::split(sequence, mode);
        i = Palantir::CompoundCodon::to_digit(s, g);
    }
    else {
        stop("Unkown sequence type");
    }

    List seq = List::create(
        _["sequence"] = s,
        _["index"] = i,
        _["type"] = type,
        _["length"] = i.size()
    );
    seq.attr("class") = "Sequence";
    return seq;
}

//[[Rcpp::export]]
List sample_sequence(List model, unsigned long long length)
{
    if(!has_class(model, "SubstitutionModel")) {
        stop("Argument `model` should be of class `SubstitutionModel`");
    }

    Palantir::GeneticCode g(get_genetic_code_name());
    uvec i = Palantir::sample_sequence(model["equilibrium"], length);
    string type = model["type"];
    vector<string> s;

    if(type == "nucleotide") {
        s = Palantir::Nucleotide::to_string(i);
    }
    else if(type == "amino_acid") {
        s = Palantir::AminoAcid::to_string(i);
    }
    else if(type == "codon") {
        s = Palantir::Codon::to_string(i, g);
    }
    else if(type == "codon_pair") {
        s = Palantir::CodonPair::to_string(i, g);
    }
    else if(type == "compound_codon") {
        s = Palantir::CompoundCodon::to_string(i, g);
    }
    else {
        stop("Unkown model type");
    }

    List seq = List::create(
        _["sequence"] = s,
        _["index"] = i,
        _["type"] = type,
        _["length"] = i.size()
    );
    seq.attr("class") = "Sequence";
    return seq;
}
