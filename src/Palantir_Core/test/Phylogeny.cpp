#include "test.hpp"

TEST_CASE("Phylogeny") {
    string newick(test_tree);
    Phylogeny tree(newick);

    SECTION("Iterating from root") {

        vector<reference_wrapper<const Phylogeny::Node>> nodes = tree.traversal();
        REQUIRE(nodes.size() == tree.n_nodes);

        vector<ullong> indeces;
        for(const Phylogeny::Node& node : nodes) {
            indeces.push_back(node.index);
        }
        REQUIRE(is_sorted(indeces.begin(), indeces.end()));
    }

    Phylogeny simple(test_simple_tree);
    cout << simple.to_string() << endl;
    cout << simple.to_JSON() << endl;

    vector<reference_wrapper<const Phylogeny::Node> > leaves = simple.leaf_traversal();
    vector<string> labels(leaves.size());
    for(ullong i = 0; i < leaves.size(); i++) {
        const Phylogeny::Node& leaf = leaves[i].get();
        labels[i] = leaf.label;
    }
    vector<string> l({"one", "two", "three", "four"});
    REQUIRE(labels == l);
}
