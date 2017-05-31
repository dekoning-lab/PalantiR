#include "GeneticCode.hpp"
using namespace Palantir;

/* Nucleotide */ 
const vector<string> Nucleotide::Nucleotides = chars("TCAG");
const ullong Nucleotide::size = 4;

string Nucleotide::to_string(ullong state)
{
    if(state >= Nucleotides.size()) {
        throw logic_error("Unexpected integer in nucleotide sequence " + ::to_string(state));
    }
    return Nucleotides[state];
}

vector<string> Nucleotide::to_string(const uvec& states)
{
    vector<string> nucleotides(states.size());
    for(ullong i = 0; i < states.size(); i++) {
        nucleotides[i] = to_string(states[i]);
    }
    return nucleotides;
}

ullong Nucleotide::to_digit(const string& state)
{
    auto index = ::find(Nucleotides.begin(), Nucleotides.end(), state);
    if(index == Nucleotides.end() || state.length() != 1) {
        throw logic_error("Unexpected letter in nucleotide sequence " + state);
    }
    return index - Nucleotides.begin();
}

uvec Nucleotide::to_digit(const vector<string>& states)
{
    uvec nucleotides(states.size());
    for(ullong i = 0; i < states.size(); i++) {
        nucleotides[i] = to_digit(states[i]);
    }
    return nucleotides;
}

vector<string> Nucleotide::split(const string& sequence)
{
    return chars(sequence);
}

bool Nucleotide::transition(ullong a, ullong b)
{
    return  (a == 0 && b== 1) || (a == 1 && b == 0) || (a == 2 && b == 3) || (a == 3 && b == 2);
}

uvec Nucleotide::transition(const uvec& a, const uvec& b)
{
    uvec tr(a.size());
    for(ullong i = 0; i < a.size(); i++) {
        tr[i] = transition(a[i], b[i]);
    }
    return tr;
}

/* Amino Acid */

const vector<string> AminoAcid::AminoAcids = chars("ACDEFGHIKLMNPQRSTVWY");
const ullong AminoAcid::size = 20;

string AminoAcid::to_string(ullong state)
{
    if(state >= AminoAcids.size()) {
        throw logic_error("Unexpected integer in amino acid sequence " + ::to_string(state));
    }
    return AminoAcids[state];
}

vector<string> AminoAcid::to_string(const uvec& states)
{
    vector<string> amino_acids(states.size());
    for(ullong i = 0; i < states.size(); i++) {
        amino_acids[i] = to_string(states[i]);
    }
    return amino_acids;
}

ullong AminoAcid::to_digit(const string& state)
{
    auto index = ::find(AminoAcids.begin(), AminoAcids.end(), state);
    if(index == AminoAcids.end() || state.length() != 1) {
        throw logic_error("Unexpected letter in amino acid sequence " + state);
    }
    return index - AminoAcids.begin();
}

uvec AminoAcid::to_digit(const vector<string>& states)
{
    uvec amino_acids(states.size());
    for(ullong i = 0; i < states.size(); i++) {
        amino_acids[i] = to_digit(states[i]);
    }
    return amino_acids;
}

vector<string> AminoAcid::split(const string& sequence)
{
    return chars(sequence);
}

/* Codon */

Codon::Codon(ullong i, ullong n_1, ullong n_2, ullong n_3, ullong aa):
    index(i),
    nucleotides({n_1, n_2, n_3}),
    amino_acid(aa) {}

string Codon::sequence() const
{
    stringstream b;
    for(ullong i = 0; i < 3; i++) {
        b << Nucleotide::Nucleotides[nucleotides[i]];
    }
    return b.str();
}

string Codon::to_string(ullong state, const GeneticCode& g)
{
    if (state >= g.size) {
        throw logic_error("Unexpected integer in codon sequence " + ::to_string(state));
    }
    return g[state].sequence();
}

vector<string> Codon::to_string(const uvec& states, const GeneticCode& g)
{
    vector<string> codons(states.size());
    for(ullong i = 0; i < states.size(); i++) {
        codons[i] = to_string(states[i], g);
    }
    return codons;
}

ullong Codon::to_digit(const string& state, const GeneticCode& g)
{
    const vector<string>& sequence = g.sequence;
    auto index = find(sequence.begin(), sequence.end(), state);
    if (index == sequence.end()) {
        throw logic_error("Unexpected string in codon sequence " + state);
    }
    return index - sequence.begin();
}

uvec Codon::to_digit(const vector<string>& states, const GeneticCode& g)
{
    uvec codons(states.size());
    for(ullong i = 0; i < states.size(); i++) {
        codons[i] = to_digit(states[i], g);
    }
    return codons;
}

vector<string> Codon::split(const string& sequence)
{
    if (sequence.length() % 3 != 0) {
        throw logic_error("Codon sequence out-of-frame");
    }
    size_t n = sequence.length() / 3;
    vector<string> strings(n);
    for(ullong i = 0; i < n; i++) {
        strings[i] = sequence.substr(i * 3, 3);
    }
    return strings;
}

string Codon::to_amino_acid(ullong state, const GeneticCode& g)
{
    return AminoAcid::AminoAcids[g[state].amino_acid];
}

vector<string> Codon::to_amino_acid(const uvec& states, const GeneticCode& g)
{
    vector<string> aa(states.size());
    for(ullong i = 0; i < states.size(); i++) {
        aa[i] = to_amino_acid(states[i], g);
    }
    return aa;
}

ullong Codon::_distance(const Codon& a, const Codon& b) 
{
    ullong d = 0;
    for (ullong i = 0; i < 3; i++) {
        if (a.nucleotides[i] != b.nucleotides[i]) d++;
    }
    return d;
}

ullong Codon::distance(ullong a, ullong b, const GeneticCode& g)
{
    return _distance(g[a], g[b]);
}

uvec Codon::distance(const uvec& a, const uvec& b, const GeneticCode& g)
{
    uvec d(a.size());
    for(ullong i = 0; i < a.size(); i++) {
        d[i] = distance(a[i], b[i], g);
    }
    return d;
}

// Requires distance(a, b) >= 1
pair<ullong, ullong> Codon::_substitution(const Codon& a, const Codon& b)
{
    for(ullong i = 0; i < 3; i++) {
        if(a.nucleotides[i] != b.nucleotides[i]) {
            return make_pair(a.nucleotides[i], b.nucleotides[i]);
        }
    }
    throw logic_error("No substitution between two of the same codon");
}

bool Codon::_synonymous(const Codon& a, const Codon& b)
{
    return a.amino_acid == b.amino_acid;
}

bool Codon::synonymous(ullong a, ullong b, const GeneticCode& g)
{
    return _synonymous(g[a], g[b]);
}

uvec Codon::synonymous(const uvec& a, const uvec& b, const GeneticCode& g)
{
    uvec s(a.size());
    for(ullong i = 0; i < a.size(); i++) {
        s[i] = synonymous(a[i], b[i], g);
    }
    return s;
}

/* Genetic Code */
const map<string,string> GeneticCode::GeneticCode_table = { // {{{
    // The Standard Code (transl_table=1)
    {"Standard nuclear", "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"},
    // The Vertebrate Mitochondrial Code (transl_table=2)
    {"Vertebrate mitochondrial", "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG"},
    // The Yeast Mitochondrial Code (transl_table=3)
    {"Yeast mitochondrial", "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG"},
    // The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code (transl_table=4)
    {"Protozoan mitochondrial", "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"},
    // The Invertebrate Mitochondrial Code (transl_table=5)
    {"Invertebrate mitochondrial", "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG"},
    // The Ciliate, Dasycladacean and Hexamita Nuclear Code (transl_table=6)
    {"Ciliate nuclear", "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"},
    // The Echinoderm and Flatworm Mitochondrial Code (transl_table=9)
    {"Echinoderm mitochondrial", "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG"},
    // The Euplotid Nuclear Code (transl_table=10)
    {"Euplotid nuclear", "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"},
    // The Bacterial, Archaeal and Plant Plastid Code (transl_table=11)
    {"Plastid", "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"},
    // The Alternative Yeast Nuclear Code (transl_table=12)
    {"Alternative yeast nuclear", "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"},
    // The Ascidian Mitochondrial Code (transl_table=13)
    {"Ascidian mitochondrial", "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG"},
    // The Alternative Flatworm Mitochondrial Code (transl_table=14)
    {"Flatworm mitochondrial", "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG"},
    // The Chlorophycean Mitochondrial Code (transl_table=16)
    {"Chlorophycean mitochondrial", "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"},
    // The Trematode Mitochondrial Code (transl_table=21)
    {"Trematode mitochondrial", "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG"},
    // The Scenedesmus obliquus Mitochondrial Code (transl_table=22)
    {"Scenedesmus mitochondrial", "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"},
    // The Thraustochytrium Mitochondrial Code (transl_table=23)
    {"Thraustochytrium mitochondrial", "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"},
    // The Pterobranchia Mitochondrial Code (transl_table=24)
    {"Pterobranchia mitochondrial", "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG"},
    // The Candidate Division SR1 and Gracilibacteria Code (transl_table=25)
    {"Gracilibacterial", "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"}};
//}}}

GeneticCode::GeneticCode(string type): size(0), type(type)
{
    auto loc = GeneticCode_table.find(type);
    if (loc == GeneticCode_table.end()) {
        throw invalid_argument("Unrecognized genetic code type in GeneticCode constructor");
    } else {
        code = loc->second;
    }

    codons = vector<Codon>();
    sequence = vector<string>();

    ullong l = 0;
    ullong m = 0;
    for (ullong i = 0; i < 4; i++) {
        for (ullong j = 0; j < 4; j++) {
            for (ullong k = 0; k < 4; k++) {
                if (code[l] != '*') {
                    string aa_s = code.substr(l, 1);
                    Codon c = Codon(m, i, j, k, AminoAcid::to_digit(aa_s));
                    codons.push_back(c);
                    sequence.push_back(c.sequence());
                    m++;
                }
                l++;
            }
        }
    }
    size = m;
}

/* Codon Pair */

string CodonPair::to_string(ullong state, const GeneticCode& g)
{
    ullong index_2 = state % g.size;
    ullong index_1 = (state - index_2) / g.size;
    return Codon::to_string(index_1, g) + "," + Codon::to_string(index_2, g);
}

vector<string> CodonPair::to_string(const uvec& states, const GeneticCode& g)
{
    vector<string> codon_pairs(states.size());
    for(ullong i = 0; i < states.size(); i++) {
        codon_pairs[i] = to_string(states[i], g);
    }
    return codon_pairs;
}

ullong CodonPair::to_digit(const string& state, const GeneticCode& g)
{
    string state_1 = state.substr(0, 3);
    string state_2 = state.substr(4, 3);
    return (Codon::to_digit(state_1, g) * g.size) + Codon::to_digit(state_2, g);
}

uvec CodonPair::to_digit(const vector<string>& states, const GeneticCode& g)
{
    uvec codon_pairs(states.size());
    for(ullong i = 0; i < states.size(); i++) {
        codon_pairs[i] = to_digit(states[i], g);
    }
    return codon_pairs;
}

ullong CodonPair::first(ullong state, const GeneticCode& g)
{
    ullong index_2 = state % g.size;
    return (state - index_2) / g.size;
}

uvec CodonPair::first(const uvec& states, const GeneticCode& g)
{
    uvec codon_pairs(states.size());
    for(ullong i = 0; i < states.size(); i++) {
        codon_pairs[i] = first(states[i], g);
    }
    return codon_pairs;
}

ullong CodonPair::second(ullong state, const GeneticCode& g)
{
    return state % g.size;
}

uvec CodonPair::second(const uvec& states, const GeneticCode& g)
{
    uvec codon_pairs(states.size());
    for(ullong i = 0; i < states.size(); i++) {
        codon_pairs[i] = second(states[i], g);
    }
    return codon_pairs;
}

vector<string> CodonPair::split(const string& sequence)
{
    if(sequence.length() % 6 != 0) {
        throw logic_error("Codon pair sequence out-of-frame");
    }
    ullong n = sequence.length() / 6;
    vector<string> strings(n);
    for(ullong i = 0; i < n; i++) {
        string a = sequence.substr(i * 6, 3);
        string b = sequence.substr((i * 6) + 3, 3);
        strings[i]  = a + "," + b;
    }
    return strings;
}

ullong CodonPair::index(const Codon& c, const Codon& d, const GeneticCode& g)
{
    return (c.index * g.size) + d.index;
}

/* Compound Codon */

string CompoundCodon::to_string(ullong state, const GeneticCode& g)
{
    ullong index = state % g.size;
    ullong mode = (state - index) / g.size;
    return Codon::to_string(index, g) + "," + ::to_string(mode);
}

vector<string> CompoundCodon::to_string(const uvec& states, const GeneticCode& g)
{
    vector<string> compound_codons(states.size());
    for(ullong i = 0; i < states.size(); i++) {
        compound_codons[i] = to_string(states[i], g);
    }
    return compound_codons;
}

ullong CompoundCodon::to_digit(const string& state, const GeneticCode& g)
{
    string codon = state.substr(0, 3);
    ullong mode = ::stoull(state.substr(4));
    return (mode * g.size) +  Codon::to_digit(codon, g);
}

uvec CompoundCodon::to_digit(const vector<string>& states, const GeneticCode& g)
{
    uvec compound_codons(states.size());
    for(ullong i = 0; i < states.size(); i++) {
        compound_codons[i] = to_digit(states[i], g);
    }
    return compound_codons;
}

vector<string> CompoundCodon::split(const string& sequence, ullong mode)
{
    vector<string> codons = Codon::split(sequence);
    vector<string> compound_codons(codons.size());

    for(ullong i = 0; i < compound_codons.size(); i++) {
        compound_codons[i] = codons[i] + "," + ::to_string(mode);
    }
    return compound_codons;
}

ullong CompoundCodon::index(const Codon& c, ullong mode, const GeneticCode& g)
{
    return (mode * g.size) + c.index;
}

ullong CompoundCodon::to_mode_index(const ullong state, const GeneticCode& g)
{
    ullong index = state % g.size;
    return (state - index) / g.size; 
}

uvec CompoundCodon::to_mode_index(const uvec& states, const GeneticCode& g)
{
    uvec compound_codons(states.size());
    for(ullong i = 0; i < states.size(); i++) {
        compound_codons[i] = to_mode_index(states[i], g);
    }
    return compound_codons;
}

ullong CompoundCodon::to_codon_index(const ullong state, const GeneticCode& g)
{
    return state % g.size;
}


uvec CompoundCodon::to_codon_index(const uvec& states, const GeneticCode& g)
{
    uvec compound_codons(states.size());
    for(ullong i = 0; i < states.size(); i++) {
        compound_codons[i] = to_codon_index(states[i], g);
    }
    return compound_codons;
}
