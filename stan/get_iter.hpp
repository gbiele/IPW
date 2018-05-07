static int iter_count = 0;

#include <iostream>
#include <fstream>
#include <string> 
using namespace std;

int get_iter(const int&id, std::ostream* pstream__) {
  
  std::string file_name = "tmp" + boost::lexical_cast<std::string>(id) + ".csv";
  //std::string file_name = "tmp" + std::to_string(id) + ".csv";
  
  int iter_count = 0;
  string line;
  
  /* Creating input filestream */ 
  ifstream file(file_name.c_str());
  while (getline(file, line))
    iter_count++;
  
  return iter_count-24;
}
