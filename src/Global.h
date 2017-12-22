#ifndef GLOBAL_H
#define GLOBAL_H

#include <string>
#include <map>
#include <vector>
#include <memory>
#include <limits>

using namespace std;

using valType_uu = std::map<unsigned,unsigned>::value_type;
using valType_us = std::map<unsigned,string>::value_type;
using valType_su = std::map<string,unsigned>::value_type;

using uint_vec_t = std::vector<unsigned int>;
using int_vec_t = std::vector<int>;
using float_vec_t = std::vector<float>;
using double_vec_t = std::vector<double>;

using uint_mat_t = std::vector< std::vector<unsigned int> >;
using int_mat_t = std::vector< std::vector<int> >;
using float_mat_t = std::vector< std::vector<float> >;
using double_mat_t = std::vector< std::vector<double> >;


#endif
