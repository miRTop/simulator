#ifndef UTILS
#define UTILS

#include <iostream>
#include <vector>
#include <string>
#include <sstream>

using namespace std;

vector<string> split(string str, char delimiter);

void reverseStr(string& str);

void path_check(std::string &path);

void delete_file(std::string file_name);
 
#endif
