#include "Phylogeny.hpp"

Palantir::Phylogeny::Phylogeny(const string &newick) : root(new Node) //Read a string and turn it into a tree
{
    string tree(newick);
    tree.erase(remove_if(tree.begin(), tree.end(), ::isspace), tree.end());
    vector<string> tokens = tokenize(tree, "(,)");
    //Initialize tree
    root->label = "root";
    root->length = 0;
    root->index = 0;

    Node *cursor = root;
    n_nodes = 1;

    //Create children
    for (auto t : tokens) {
        if (t == "(") { //New node
            Node *temp = new Node;
            temp->index = n_nodes;
            cursor->children.push_back(temp);
            temp->parent = cursor;
            temp->parent_index = cursor->index;
            cursor = temp;
            n_nodes++;

        } else if (t == ",") { //New sibling
            Node *temp = new Node;
            Node *parent = cursor->parent;
            temp->index = n_nodes;
            parent->children.push_back(temp);
            temp->parent = parent;
            temp->parent_index = parent->index;
            cursor = temp;
            n_nodes++;
        } else if (t == ")") { //End clade
            cursor = cursor->parent;
        } else { //Record node length and label
            ullong colon = t.find(":"); //find position of the colon character
            if (colon != string::npos) //if colon symbol is found
            {
                string label = t.substr(0, colon); //try to get label
                double length = stod(t.substr(colon + 1, t.length() - colon));
                cursor->length = length;
                if (label == "") //No label - give a number index
                {
                    cursor->label = ::to_string(cursor->index);
                } else //Give a label from tree
                {
                    cursor->label = label;
                }
            }
        } //
    }
}

Palantir::Phylogeny::Node::~Node() 
{
    for (Node *c : children) {
        delete c;
    }
}

Palantir::Phylogeny::~Phylogeny() 
{
    delete root;
}

ostream &Palantir::Phylogeny::Node::print(ostream &os, string prefix, bool isLast) const 
{
    os << prefix << (isLast ? "└──" : "├──") << "[" << index << "] \"" << label << "\"" << ", " << length << endl;
    if (!children.empty()) { //if node is not terminal
        for (auto c : children) {
            c->print(os, prefix + (isLast ? "   " : "│  "),
                     c == children.back()); // space if terminal, bar '|' if it has siblings
        }
    }
    return os;
}

const vector<string> Palantir::Phylogeny::tokenize(const string &str, const char *delims) const 
{
    vector<string> tokens;
    ullong start = 0;
    for (ullong end = 0; end != string::npos; end = str.find_first_of(delims, start)) {
        if (end != start) {
            tokens.push_back(str.substr(start, end - start));
        }
        tokens.push_back(str.substr(end, 1));
        start = end + 1;
    }
    return tokens;
}

void Palantir::Phylogeny::Node::traverse(vector<reference_wrapper<const Node>>& walk) const
{
    walk.push_back(*this);
    for (auto child : children) {
        child->traverse(walk);
    }
}

void Palantir::Phylogeny::Node::traverse_leaves(vector<reference_wrapper<const Node>>& walk) const
{
    if(this->is_leaf()) {
        walk.push_back(*this);
    }
    for (auto child : children) {
        child->traverse_leaves(walk);
    }
}

vector<reference_wrapper<const Palantir::Phylogeny::Node>> Palantir::Phylogeny::traversal() const
{
    vector<reference_wrapper<const Phylogeny::Node>> walk;
    root->traverse(walk);
    return (walk);
}

vector<reference_wrapper<const Palantir::Phylogeny::Node>> Palantir::Phylogeny::leaf_traversal() const
{
    vector<reference_wrapper<const Phylogeny::Node>> walk;
    root->traverse_leaves(walk);
    return (walk);
}

uvec Palantir::Phylogeny::leaf_indeces() const
{
    vector<reference_wrapper<const Phylogeny::Node> > leaves = leaf_traversal();
    uvec indeces(leaves.size());
    for(ullong i = 0; i < leaves.size(); i++) {
        const Phylogeny::Node& node = leaves[i].get();
        indeces[i] = node.index;
    }
    return indeces;
}

void Palantir::Phylogeny::Node::to_JSON(stringstream& b) const 
{
    b << "{";
    b << "\"index\":" << index << ",";
    b << "\"label\":\"" << label << "\",";
    b << "\"length\":" << length;

    if (!is_leaf()) {
        b << ",\"children\":[";
    }

    for (auto child : children) {
        child->to_JSON(b);
        if (child != children.back()) {
            b << ",";
        } else {
            b << "]";
        }
    }
    b << "}";
}

string Palantir::Phylogeny::to_JSON() const 
{
    stringstream b;
    root->to_JSON(b);
    return b.str();
}

vector<Palantir::IntervalHistory> Palantir::Phylogeny::to_intervals(const Palantir::Phylogeny& mode_tree) const 
{
    vector<reference_wrapper<const Phylogeny::Node> > mode_nodes = mode_tree.traversal();
    vector<reference_wrapper<const Phylogeny::Node> > nodes = traversal();

    vector<IntervalHistory> tree_intervals(n_nodes);
    for(const Phylogeny::Node& node : nodes) {
        ullong i = node.index;
        if(!node.is_root()) {
            const Phylogeny::Node& mode_node = mode_nodes[i].get();
            ullong state = (ullong)mode_node.length;
            tree_intervals[i].append(state, 0, node.length);
        }
    }

    return tree_intervals;
}
