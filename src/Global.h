#ifndef GLOBAL_H
#define GLOBAL_H

#include <string>
#include <map>
#include <vector>
#include <memory>
#include <limits>

using namespace std;

using uu_map_t = std::map<unsigned int,unsigned int>;
using us_map_t = std::map<unsigned int,string>;
using su_map_t = std::map<string,unsigned int>;

using uint_vec_t = std::vector<unsigned int>;
using int_vec_t = std::vector<int>;
using float_vec_t = std::vector<float>;
using double_vec_t = std::vector<double>;

using uint_mat_t = std::vector< std::vector<unsigned int> >;
using int_mat_t = std::vector< std::vector<int> >;
using float_mat_t = std::vector< std::vector<float> >;
using double_mat_t = std::vector< std::vector<double> >;


#endif
