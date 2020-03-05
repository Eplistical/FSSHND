#include "misc/ioer.hpp"
#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <utility>
using namespace std;

int main(int argc, char** argv) {
    vector<pair<double, string>> v{ make_pair(1.1, "haha"), make_pair(2.2, "hehe")};

    vector<
        pair<
            pair<vector<string>, double>
            , 
            pair<vector<double>, string>
            >
        > 
        vv(3, 
            make_pair( 
                make_pair(
                    vector<string> { "a", "b" }, 3.0
                    ), 
                make_pair(
                    vector<double> { 1.1, 2.2 }, "heihei"
                    )
                )
            );

    map<int, pair<double, string>> 
        mm {
            { 3, make_pair(1.1, "heihei") },
            { 2, make_pair(2.2, "hehe") },
            { 1, make_pair(3.3, "haha") },
            { 0, make_pair(4.4, "hoho") },
        };

    set<pair<double, map<string, int>>> 
        ss {
            make_pair(1.1, map<string, int> { {"hehe", 1} } ),
            make_pair(2.2, map<string, int> { {"heihei", 2}, {"haha", 42} } ),
            make_pair(3.3, map<string, int> { {"hoho", 3} } ),
        };

    unordered_set<double> sss {1,2,3,4};

    vector<vector<vector<pair<double, string>>>> a(3, vector<vector<pair<double, string>>>(5, vector<pair<double, string>>(v)));
    ioer::info(a);
    ioer::tabout(a);

    ioer::info(vv);
    ioer::tabout(vv);

    ioer::info(mm);
    ioer::tabout(mm);

    ioer::info(ss);
    ioer::tabout(ss);

    ioer::info(sss);
    ioer::tabout(111, sss);
    return 0;
}
