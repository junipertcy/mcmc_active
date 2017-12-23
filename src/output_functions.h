//
// Created by Tzu-Chi Yen on 12/24/17.
//

#ifndef MCMC_ACTIVE_OUTPUT_FUNCTIONS_H
#define MCMC_ACTIVE_OUTPUT_FUNCTIONS_H

#include <iostream>
#include "Global.h"

template<typename T>
void output_mat(T mat, std::ostream & stream=std::clog)
{
    for (unsigned int r = 0; r < mat.size(); ++ r)
    {
        for (unsigned int s = 0; s < mat[r].size(); ++s)
        {
            stream << mat[r][s] << " ";
        }
        stream << "\n";
    }
}
template<typename T>
void output_vec(T vec, std::ostream & stream=std::clog)
{
    for (auto it = vec.begin(); it != vec.end(); ++it)
    {
        stream << *it << " ";
    }
    stream << "\n";
}


#endif //MCMC_ACTIVE_OUTPUT_FUNCTIONS_H
