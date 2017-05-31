#include "Util.hpp"

void removeWhitespace(string& str)
{
    str.erase(
            remove_if(str.begin(), str.end(), ::isspace),
            str.end());
}

Match find_first_in(string& str, vector<string> patterns, ullong pos)
{
    vector<ullong> found;
    for(auto pattern : patterns) {
        found.push_back(str.find(pattern, pos));
    }
    long first_match = distance(found.begin(), min_element(found.begin(), found.end()));

    if(found[first_match] == string::npos) {
        return noMatch;
    }
    else {
        return Match(found[first_match], patterns[first_match].size());
    }

}

vector<string> tokenize(string& str, vector<string> delimiters, bool keepDelimiters, bool keepEmpty)
{
    vector<string> tokens;
    ullong start_pos = 0;
    if(str.empty()) {
        return tokens;
    }
    if(delimiters.empty()) {
        return tokens;
    }

    Match match = find_first_in(str, delimiters, start_pos);

    while(match != noMatch) {

        tokens.push_back(str.substr(start_pos, match.begin - start_pos));
        if(keepDelimiters) {
            tokens.push_back(str.substr(match.begin, match.length));
        }

        start_pos = match.begin + match.length;
        match = find_first_in(str, delimiters, start_pos);
    }
    // Push the rest
    tokens.push_back(str.substr(start_pos));

    if(!keepEmpty) {
        tokens.erase(
                remove_if(tokens.begin(), tokens.end(), [](const string& s)->bool { return s.empty(); }),
                tokens.end());
    }

    return tokens;
}

string file_to_string(const string path) {
    ifstream t(path);
    if (!t.is_open()) {
        throw invalid_argument("Can not open file '" + path + "'");
    }
    stringstream buffer;
    buffer << t.rdbuf();
    return buffer.str();
}

vector<string> chars(const string& s) 
{
    vector<string> ch;
    for(const char& c : s) {
        ch.push_back(string(1, c));
    }
    return ch;
}
