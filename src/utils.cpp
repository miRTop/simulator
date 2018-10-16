#include <iostream>
#include <vector>
#include <string>
#include <sstream>

using namespace std;

//Function to split a string by specified delimeter.
vector<string> split(string str, char delimiter) {
  vector<string> internal;
  stringstream ss(str); // Turn the string into a stream.
  string tok;
 
  while(getline(ss, tok, delimiter)) {
    internal.push_back(tok);
  }
 
  return internal;
}

//Function to reverse a string.
void reverseStr(string& str) 
{ 
    int n = str.length(); 
  
    // Swap character starting from two 
    // corners 
    for (int i = 0; i < n / 2; i++) 
        swap(str[i], str[n - i - 1]); 
} 

void path_check(std::string &path){ //FIXME check if path exists 
  
  std::string last = path.substr(path.size(),1);
  if(last != "/"){
    path = path + "/";
  }
  
}

void delete_file(std::string file_name){
  
 if(remove(file_name.c_str()) != 0){
    perror("Error deleting a file\n");
 }
 
}
