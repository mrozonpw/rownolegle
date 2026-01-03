#ifndef ALGORITHM_H
#define ALGORITHM_H
#pragma once
#include <vector>
#include <utility>
#include <cstdint>
#include "functions.h"

// sekwencyjne SA
std::pair<std::vector<double>, double>
perform_sequential_algorithm(const calc_function_t& calc_value,
                             std::vector<double> x0,
                             uint32_t n,
                             int a, int b);

// wielowątkowe SA
std::pair<std::vector<double>, double>
perform_threaded_algorithm(const calc_function_t& calc_value,
                           uint32_t n,
                           int a, int b,
                           int num_threads);


#endif //ALGORITHM_H
