#ifndef GeneticCode_hpp
#define GeneticCode_hpp

#include "Util.hpp"

using namespace std;

namespace Palantir 
{

    class Nucleotide;
    class AminoAcid;
    class Codon;
    class GeneticCode;
    class CodonPair;
    class CompoundCodon;

    class Nucleotide
    {
        public:
            const static vector<string> Nucleotides;
            const static ullong size;

            static string to_string(ullong state);
            static ullong to_digit(const string& state);

            static vector<string> to_string(const uvec& states);
            static uvec to_digit(const vector<string>& states);

            static vector<string> split(const string& sequence);

            static bool transition(ullong a, ullong b);
            static uvec transition(const uvec& a, const uvec& b);
    };

    class AminoAcid
    {
        public:
            const static vector<string> AminoAcids;
            const static ullong size;

            static string to_string(ullong state);
            static ullong to_digit(const string& state);

            static vector<string> split(const string& sequence);

            static vector<string> to_string(const uvec& states);
            static uvec to_digit(const vector<string>& states);
    };

    class Codon
    {
        public:
            ullong index;
            vector<ullong> nucleotides;
            ullong amino_acid;

            static string to_string(ullong state, const GeneticCode& g);
            static ullong to_digit(const string& state, const GeneticCode& g);

            static vector<string> to_string(const uvec& states, const GeneticCode& g);
            static uvec to_digit(const vector<string>& states, const GeneticCode& g);

            static vector<string> split(const string& sequence);

            static string to_amino_acid(ullong state, const GeneticCode& g);
            static vector<string> to_amino_acid(const uvec& states, const GeneticCode& g);


            Codon(ullong i, ullong n_1, ullong n_2, ullong n_3, ullong aa);
            string sequence() const;

            friend inline bool operator==(const Codon& lhs, const Codon& rhs) {
                return lhs.index == rhs.index;
            }
            friend inline bool operator!=(const Codon& lhs, const Codon& rhs) {
                return !(lhs == rhs);
            }

            // with codon instance
            static pair<ullong, ullong> _substitution(const Codon& a, const Codon& b);
            static ullong _distance(const Codon& a, const Codon& b);
            static bool _synonymous(const Codon& a, const Codon& b);

            // with digits
            static ullong distance(ullong a, ullong b, const GeneticCode& g);
            static uvec distance(const uvec& a, const uvec& b, const GeneticCode& g);

            static bool synonymous(ullong a, ullong b, const GeneticCode& g);
            static uvec synonymous(const uvec& a, const uvec& b, const GeneticCode& g);

    };

    class GeneticCode
    {
        public:
            const static map<string, string> GeneticCode_table;

            const Codon& operator[](ullong i) const { return codons[i]; }
            ullong size;

            vector<Codon> codons;
            vector<string> sequence;

            GeneticCode(string type);

            vector<Codon>::const_iterator begin() const { return codons.begin(); }
            vector<Codon>::const_iterator end() const { return codons.end(); }

        private:
            string code;
            string type;
    };

    class CodonPair
    {
        public:
            static string to_string(ullong state, const GeneticCode& g);
            static ullong to_digit(const string& state, const GeneticCode& g);

            static vector<string> to_string(const uvec& states, const GeneticCode& g);
            static uvec to_digit(const vector<string>& states, const GeneticCode& g);

            static ullong first(ullong state, const GeneticCode& g);
            static uvec first(const uvec& states, const GeneticCode& g);

            static ullong second(ullong state, const GeneticCode& g);
            static uvec second(const uvec& states, const GeneticCode& g);

            static vector<string> split(const string& sequence);

            static ullong index(const Codon& c, const Codon& d, const GeneticCode& g);
    };

    class CompoundCodon
    {
        public:
            static string to_string(ullong state, const GeneticCode& g);
            static ullong to_digit(const string& state, const GeneticCode& g);

            static vector<string> to_string(const uvec& states, const GeneticCode& g);
            static uvec to_digit(const vector<string>& states, const GeneticCode& g);

            static vector<string> split(const string& sequence, ullong mode);

            static ullong index(const Codon& c, ullong mode, const GeneticCode& g);

            static ullong to_mode_index(const ullong state, const GeneticCode& g);
            static uvec to_mode_index(const uvec& states, const GeneticCode& g);

            static ullong to_codon_index(const ullong state, const GeneticCode& g);
            static uvec to_codon_index(const uvec& states, const GeneticCode& g);
    };
}

#endif /* GeneticCode_hpp */
