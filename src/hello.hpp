#ifndef __HELLO_HPP__
#define __HELLO_HPP__

#warning "bla blubb"

#include <boost/foreach.hpp>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#define foreach BOOST_FOREACH

using namespace std;
struct World
{
	void geti(int& i){ i = 2; }
    void set(std::string msg) { ss.push_back(msg); }
    std::string greetx(World& x) {
		stringstream sstream;
		foreach( std::string& s, x.ss){
			sstream<<s<<endl;
		}
		return sstream.str();
	}
    std::string greet() {
		stringstream sstream;
		foreach( std::string& s, ss){
			sstream<<s<<endl;
		}
		return sstream.str();
	}
	std::vector<string> ss;
	virtual string f(){
		return "in parent";
	}
};

struct W2 : public World{
	virtual string f(){
		return "in child";
	}
};

#endif
