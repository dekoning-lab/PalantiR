#ifndef Phylogeny_hpp
#define Phylogeny_hpp

#include "Util.hpp"
#include "IntervalHistory.hpp"

namespace Palantir {
    class Phylogeny {
        public:
            Phylogeny(const string &newick);

            ~Phylogeny();

            friend ostream &operator<<(ostream &os, const Phylogeny &p) {
                p.root->print(os, "", true);
                return os;
            }

            string to_string() {
                ostringstream b;
                root->print(b, "", true);
                return b.str();
            }

            class Node {
                private:
                    Node *parent;
                    vector<Node *> children;
                    void traverse(vector<reference_wrapper<const Node>>& walk) const;
                    void traverse_leaves(vector<reference_wrapper<const Node>>& walk) const;
                    void to_JSON(stringstream& b) const;

                public:
                    friend class Phylogeny;
                    double length;
                    ullong index;
                    ullong parent_index;
                    string label;

                    ostream &print(ostream &os, string prefix, bool isLast) const;

                    Node() : parent(nullptr), children(vector<Node *>()), length(0), index(0), parent_index(0), label("") {}

                    bool is_leaf() const { return children.empty(); }
                    bool is_root() const { return parent == nullptr; }

                    ~Node();
            };

            ullong n_nodes;

            vector<reference_wrapper<const Node>> traversal() const;
            vector<reference_wrapper<const Node>> leaf_traversal() const;
            uvec leaf_indeces() const;

            string to_JSON() const;

            // Use branch length as interval indicators
            vector<IntervalHistory> to_intervals(const Phylogeny& mode_tree) const;

        private:
            const vector<string> tokenize(const string &newick, const char *delim) const;
            Node *root;
    };
}

#endif /* Phylogeny_hpp */
