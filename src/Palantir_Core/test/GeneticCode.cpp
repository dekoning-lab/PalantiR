#include "test.hpp"

TEST_CASE("GeneticCode")
{
    GeneticCode g("Standard nuclear");

    REQUIRE(g.size == 61);
    REQUIRE(g.codons.size() == 61);
    REQUIRE(g.sequence.size() == 61);

    for(auto c : g) {
        REQUIRE(c.sequence().size() == 3);
    }

    GeneticCode mito("Vertebrate mitochondrial");
    REQUIRE(mito.size == 60);
}

TEST_CASE("Nucleotide")
{

    vector<string> nuc = Palantir::Nucleotide::split("CCTAGT");
    uvec nuc_i = Palantir::Nucleotide::to_digit(nuc);
    vector<string> nuc_s = Palantir::Nucleotide::to_string(nuc_i);
    REQUIRE(nuc_s == nuc);

    REQUIRE_THROWS_WITH(Palantir::Nucleotide::to_digit("H"), Catch::Equals("Unexpected letter in nucleotide sequence H"));
    REQUIRE_THROWS_WITH(Palantir::Nucleotide::to_string(5), Catch::Equals("Unexpected integer in nucleotide sequence 5"));
}

TEST_CASE("AminoAcid")
{
    vector<string> aa = Palantir::AminoAcid::split("AKTWYCGRSMW");
    uvec aa_i = Palantir::AminoAcid::to_digit(aa);
    vector<string> aa_s = Palantir::AminoAcid::to_string(aa_i);
    REQUIRE(aa_s == aa);
    REQUIRE_THROWS_WITH(Palantir::AminoAcid::to_digit("O"), Catch::Equals("Unexpected letter in amino acid sequence O"));
    REQUIRE_THROWS_WITH(Palantir::AminoAcid::to_string(20), Catch::Equals("Unexpected integer in amino acid sequence 20"));
}

TEST_CASE("Codon")
{
    GeneticCode g("Standard nuclear");
    vector<string> codon = Palantir::Codon::split("ACTCATCGGGAGACTTATCGAAAACATAGC");
    uvec codon_i = Palantir::Codon::to_digit(codon, g);
    vector<string> codon_s = Palantir::Codon::to_string(codon_i, g);
    REQUIRE(codon_s == codon);
    REQUIRE_THROWS_WITH(Palantir::Codon::to_digit("TGA", g), Catch::Equals("Unexpected string in codon sequence TGA"));
    REQUIRE_THROWS_WITH(Palantir::Codon::to_string(64, g), Catch::Equals("Unexpected integer in codon sequence 64"));
    REQUIRE_THROWS_WITH(Palantir::Codon::split("ACTC"), Catch::Equals("Codon sequence out-of-frame"));
}

TEST_CASE("CodonPair")
{
    GeneticCode g("Standard nuclear");
    vector<string> pairs = Palantir::CodonPair::split("ACTACTCATACTCGGGGCGAGCACACTTACTATTATCGAACGAAAAAACATCATAGCACG");
    vector<string> codon_pairs = {"ACT,ACT", "CAT,ACT", "CGG,GGC", "GAG,CAC", "ACT,TAC", "TAT,TAT", "CGA,ACG", "AAA,AAA", "CAT,CAT", "AGC,ACG"};
    REQUIRE(pairs == codon_pairs);
    uvec codon_pairs_i = Palantir::CodonPair::to_digit(codon_pairs, g);
    codon_pairs_i.print();
    vector<string> codon_pairs_s = Palantir::CodonPair::to_string(codon_pairs_i, g);
    REQUIRE_THROWS_WITH(Palantir::CodonPair::to_digit("TAG,TTT", g), Catch::Equals("Unexpected string in codon sequence TAG"));
    REQUIRE(codon_pairs_s == codon_pairs);

}

TEST_CASE("CompoundCodon")
{
    GeneticCode g("Standard nuclear");
    vector<string> compound_states = Palantir::CompoundCodon::split("ACTCATCGGGAGACTTATCGAAAACATAGC", 2);
    vector<string> compound_codons = {"ACT,2", "CAT,2", "CGG,2", "GAG,2", "ACT,2", "TAT,2", "CGA,2", "AAA,2", "CAT,2", "AGC,2"};
    REQUIRE(compound_states == compound_codons);
    uvec compound_states_i = Palantir::CompoundCodon::to_digit(compound_states, g);
    compound_states_i.print();
    vector<string> compound_states_s = Palantir::CompoundCodon::to_string(compound_states_i, g);
    REQUIRE_THROWS_WITH(Palantir::CompoundCodon::to_digit("TAG,3", g), Catch::Equals("Unexpected string in codon sequence TAG"));
    REQUIRE(compound_states_s == compound_states);

    vector<string> z = Palantir::CompoundCodon::split("TTC", 234);
    uvec i = Palantir::CompoundCodon::to_digit(z, g);
    REQUIRE(i[0] == 14275); 
    REQUIRE(Palantir::CompoundCodon::to_codon_index(i[0], g) == 1);
    REQUIRE(Palantir::CompoundCodon::to_mode_index(i[0], g) == 234);
}
