#include "test.hpp"

TEST_CASE("Removing characters", "stringUtil")
{
    string example("hello world: the testening");

    SECTION("removing whitespace")
    {
        string t(example);
        removeWhitespace(t);
        REQUIRE(t == "helloworld:thetestening");
    }
}

TEST_CASE("find_first_in")
{
    string example("hello world: the testening");

    SECTION("Matches are comparable")
    {
        Match one(12, 14);
        Match two(12, 14);
        Match three(12, 13);
        REQUIRE(one == two);
        REQUIRE(one != three);
    }

    SECTION("Matches are printable")
    {
        Match one(12, 2);
        stringstream s;
        s << one;
        REQUIRE(s.str() == "Match(begin=12, length=2)");
    }

    SECTION("finding first delimiter")
    {
        vector<string> delims({"world", "what"});
        Match match(6, 5);
        auto found = find_first_in(example, delims);
        REQUIRE(found == match);
    }

    SECTION("finding first occurrence")
    {
        vector<string> delimiters({"or", "el"});
        Match match(1, 2);
        REQUIRE(find_first_in(example, delimiters) == match);
    }

    SECTION("finding second delim")
    {
        vector<string> delims({"what", "the"});
        Match match(13, 3);
        REQUIRE(find_first_in(example, delims) == match);
    }

    SECTION("finding no delims")
    {
        vector<string> delims({"what", "banana"});
        REQUIRE(find_first_in(example, delims) == noMatch);
    }

    SECTION("finding no delims - late start")
    {
        vector<string> delims({"what", "hello"});
        REQUIRE(find_first_in(example, delims, 1) == noMatch);
    }
}

TEST_CASE("tokenize")
{
    string example("token:token,token:token");
    SECTION("Simple tokenization")
    {
        vector<string> tokens = tokenize(example, {":", ","}, false);
        REQUIRE(tokens[0] == "token");
        REQUIRE(tokens[1] == "token");
        REQUIRE(tokens[2] == "token");
        REQUIRE(tokens[3] == "token");
        REQUIRE(tokens.size() == 4);
    }
    SECTION("Tokenization with delimeters")
    {
        vector<string> tokens = tokenize(example, {":", ","});
        REQUIRE(tokens[0] == "token");
        REQUIRE(tokens[1] == ":");
        REQUIRE(tokens[2] == "token");
        REQUIRE(tokens[3] == ",");
        REQUIRE(tokens[4] == "token");
        REQUIRE(tokens[5] == ":");
        REQUIRE(tokens[6] == "token");
        REQUIRE(tokens.size() == 7);
    }
    SECTION("Tokenization with strings")
    {
        vector<string> tokens = tokenize(example, {"token"}, false);
        REQUIRE(tokens[0] == "");
        REQUIRE(tokens[1] == ":");
        REQUIRE(tokens[2] == ",");
        REQUIRE(tokens[3] == ":");
        REQUIRE(tokens[4] == "");
        REQUIRE(tokens.size() == 5);
    }
    SECTION("Tokenization with strings - no empties")
    {
        vector<string> tokens = tokenize(example, {"token"}, false, false);
        REQUIRE(tokens[0] == ":");
        REQUIRE(tokens[1] == ",");
        REQUIRE(tokens[2] == ":");
        REQUIRE(tokens.size() == 3);
    }
}
