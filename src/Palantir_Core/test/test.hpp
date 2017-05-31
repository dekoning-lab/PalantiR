#ifndef Test_hpp
#define Test_hpp

#include <armadillo>

#include "../include/catch.hpp"

#include "../Util.hpp"
#include "../GeneticCode.hpp"

#include "../SubstitutionHistory.hpp"
#include "../IntervalHistory.hpp"
#include "../SiteSimulation.hpp"

#include "../Palantir.hpp"
#include "../HasegawaKishinoYano.hpp"
#include "../MutationSelection.hpp"
#include "../GeneralTimeReversible.hpp"
#include "../MarkovModel.hpp"
#include "../CoEvolution.hpp"
#include "../MarkovModulated.hpp"
#include "../Simulate.hpp"

using namespace Palantir;

// Fitness profile 4, N=1000
static const string test_fitness_1 = "0.9995424 0.9982223 1.0000000 1.0000000 0.9993214 0.9995571 0.9996736 0.9987027 0.9998721 0.9997290 0.9995571 1.0000000 0.9989922 0.9997263 0.9992135 0.9996857 1.0000000 0.6365381 0.9988228 0.9997598";
// Fitness profile 5, N=1000
static const string test_fitness_2 = "0.9990741 0.6355499 0.6358824 0.6359413 0.9989646 0.6365381 0.6363081 1.0000000 0.6365066 0.9999538 0.9995469 0.6362290 0.6363198 0.6362492 0.6359919 0.6362086 0.9995317 0.6365381 1.0000000 0.6365381";

static const string test_equilibrium_1 = "0.05154639 0.01030928 0.18556701 0.16494845 0.02061856 0.02061856 0.04123711 0.03092784 0.04123711 0.05154639 0.02061856 0.09278351 0.01030928 0.05154639 0.03092784 0.05154639 0.05154639 0.00000000 0.04123711 0.03092784";

static const string test_equilibrium_2 = "0.02020202 0.00000000 0.00000000 0.00000000 0.01010101 0.00000000 0.00000000 0.41414141 0.00000000 0.08080808 0.02020202 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.02020202 0.00000000 0.43434343 0.00000000";

static const string test_tree = "((((((((((((((Human:0.0067,Pan_troglodytes:0.006667):0.00225,Gorilla_gorilla:0.008825):0.00968,Pongo_abelii:0.018318):0.00717,Nomascus_leucogenys:0.025488):0.00717,(Macaca_mulatta:0.007853,Papio_hamadryas:0.007637):0.029618):0.021965,Callithrix_jacchus:0.066131):0.05759,Tarsius_syrichta:0.137823):0.011062,(Microcebus_murinus:0.092749,Otolemur_garnettii:0.129725):0.035463):0.015494,Tupaia_belangeri:0.186203):0.004937,(((((Mus_musculus:0.084509,Rattus_norvegicus:0.091589):0.197773,Dipodomys_ordii:0.211609):0.022992,Cavia_porcellus:0.225629):0.01015,Ictidomys_tridecemlineatus:0.148468):0.025746,(Oryctolagus_cuniculus:0.114227,Ochotona_princeps:0.201069):0.101463):0.015313):0.020593,((((Vicugna_pacos:0.107275,(Tursiops_truncatus:0.064688,(Bos_taurus:0.061796,Ovis_aries:0.061796):0.061796):0.025153):0.0201675,Sus_scrofa:0.079):0.0201675,((Equus_caballus:0.109397,(Felis_catus:0.098612,((Ailuropoda_melanoleuca:0.025614,Mustela_putorius_furo:0.0256):0.0256145,Canis_familiaris:0.051229):0.051229):0.049845):0.006219,(Myotis_lucifugus:0.14254,Pteropus_vampyrus:0.113399):0.033706):0.004508):0.011671,(Erinaceus_europaeus:0.221785,Sorex_araneus:0.269562):0.056393):0.021227):0.023664,(((Loxodonta_africana:0.082242,Procavia_capensis:0.155358):0.02699,Echinops_telfairi:0.245936):0.049697,(Dasypus_novemcinctus:0.116664,Choloepus_hoffmanni:0.096357):0.053145):0.006717):0.234728,(Monodelphis_domestica:0.125686,(Macropus_eugenii:0.101004,Sarcophilus_harrisii:0.101004):0.021004):0.2151):0.071664,Ornithorhynchus_anatinus:0.456592);";

static const string test_tree_mode = "((((((((((((((Human:1,Pan_troglodytes:1):1,Gorilla_gorilla:1):1,Pongo_abelii:1):1,Nomascus_leucogenys:1):1,(Macaca_mulatta:1,Papio_hamadryas:1):1):1,Callithrix_jacchus:1):1,Tarsius_syrichta:1):1,(Microcebus_murinus:1,Otolemur_garnettii:1):1):1,Tupaia_belangeri:1):1,(((((Mus_musculus:0,Rattus_norvegicus:0):0,Dipodomys_ordii:0):0,Cavia_porcellus:0):0,Ictidomys_tridecemlineatus:0):0,(Oryctolagus_cuniculus:0,Ochotona_princeps:0):0):0):0,((((Vicugna_pacos:0,(Tursiops_truncatus:0,(Bos_taurus:0,Ovis_aries:0):0):0):0,Sus_scrofa:0.):0,((Equus_caballus:0,(Felis_catus:0,((Ailuropoda_melanoleuca:0,Mustela_putorius_furo:0):0,Canis_familiaris:0):0):0):0,(Myotis_lucifugus:0,Pteropus_vampyrus:0):0):0):0,(Erinaceus_europaeus:0,Sorex_araneus:0):0):0):0,(((Loxodonta_africana:0,Procavia_capensis:0):0,Echinops_telfairi:0):0,(Dasypus_novemcinctus:0,Choloepus_hoffmanni:0):0):0):0,(Monodelphis_domestica:0,(Macropus_eugenii:0,Sarcophilus_harrisii:0):0):0):0,Ornithorhynchus_anatinus:0);";

static const string test_simple_tree = "(left(one:2,two:3):1,right(three:5,four:6):4);";

#endif //Test_hpp
