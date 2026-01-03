#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#pragma once
#include <vector>
#include <cstdint>

using calc_function_t = double (*)(const std::vector<double>&, uint32_t);

// funkcje testowe
double calc_quadratic_function(const std::vector<double>& x, uint32_t n);
double calc_woods_function(const std::vector<double>& x, uint32_t n);
double calc_powell_singular_function(const std::vector<double>& x, uint32_t n);

// punkty startowe
std::vector<double> make_quadratic_x0(uint32_t n);
std::vector<double> make_woods_x0(uint32_t n);
std::vector<double> make_powell_x0(uint32_t n);

#endif //FUNCTIONS_H
