#include <RcppArmadillo.h>

#include "Palantir_Core/GeneticCode.hpp"
#include "Palantir_Core/Palantir.hpp"

#include "RcppPalantir.hpp"

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
List as_compound(List sequence, unsigned long long mode)
{
    Palantir::GeneticCode g(get_genetic_code_name());
    if(!has_class(sequence, "Sequence")) {
        stop("Argument `sequence` should be of class `Sequence`");
    }
    vector<string> s = sequence["sequence"];
    uvec i = sequence["index"];
    string type = sequence["type"];
    unsigned long long length = sequence["length"];

    unsigned long long size;
    if(type == "nucleotide") size = Palantir::Nucleotide::size;
    else if(type == "amino_acid") size = Palantir::AminoAcid::size;
    else if(type == "codon") size = g.size;
    else if(type == "codon_pair") size = g.size * g.size;
    else if(type == "compound_codon") size = g.size;
    else stop("Unkown model type");

    i += size * mode;
    for(unsigned long long site = 0; site < length; site++) {
        s[site] += "," + ::to_string(mode);
    }

    List seq = List::create(
        _["sequence"] = s,
        _["index"] = i,
        _["type"] = "compound_codon",
        _["length"] = length
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

    if(type == "nucleotide") s = Palantir::Nucleotide::to_string(i);
    else if(type == "amino_acid") s = Palantir::AminoAcid::to_string(i);
    else if(type == "codon") s = Palantir::Codon::to_string(i, g);
    else if(type == "codon_pair") s = Palantir::CodonPair::to_string(i, g);
    else if(type == "compound_codon") s = Palantir::CompoundCodon::to_string(i, g);
    else stop("Unkown model type");

    List seq = List::create(
        _["sequence"] = s,
        _["index"] = i,
        _["type"] = type,
        _["length"] = i.size()
    );
    seq.attr("class") = "Sequence";
    return seq;
}

CharacterMatrix get_alignment(const vector<Palantir::SiteSimulation>& sims, const Palantir::Phylogeny& p, std::string type)
{
    Palantir::GeneticCode g(get_genetic_code_name());
    vector<reference_wrapper<const Palantir::Phylogeny::Node> > leaves = p.leaf_traversal();
    CharacterMatrix alignment(leaves.size(), sims.size());

    CharacterVector labels(leaves.size());
    for(unsigned long long leaf = 0; leaf < leaves.size(); leaf++) {
        const Palantir::Phylogeny::Node& node = leaves[leaf].get();
        labels[leaf] = node.label;
    }

    for(unsigned long long site = 0; site < sims.size(); site++) {
        uvec i = sims[site].leaf_states(p);
        vector<string> s;

        if(type == "nucleotide") s = Palantir::Nucleotide::to_string(i);
        else if(type == "amino_acid") s = Palantir::AminoAcid::to_string(i);
        else if(type == "codon") s = Palantir::Codon::to_string(i, g);
        else if(type == "codon_pair") s = Palantir::CodonPair::to_string(i, g);
        else if(type == "compound_codon") s = Palantir::CompoundCodon::to_string(i, g);
        else stop("Unkown model type");

        for(unsigned long long leaf = 0; leaf < leaves.size(); leaf++) {
            alignment(leaf, site) = s[leaf];
        }
    }
    rownames(alignment) = labels;
    return alignment;
}

//[[Rcpp::export]]
std::vector<std::string> as_amino_acid(std::vector<std::string>& codons)
{
    Palantir::GeneticCode g(get_genetic_code_name());
    uvec states = Palantir::Codon::to_digit(codons, g);
    return Palantir::Codon::to_amino_acid(states, g);
}
