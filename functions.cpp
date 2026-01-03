#include "functions.h"
#include <cmath>

// ===== FUNKCJE =====

double calc_quadratic_function(const std::vector<double>& x, uint32_t n)
{
    double sum = 0.0;
    for (uint32_t i = 0; i < n; ++i)
        sum += x[i] * x[i];
    return sum;
}

double calc_woods_function(const std::vector<double>& x, uint32_t n)
{
    double sum = 0.0;
    for (uint32_t i = 0; i + 3 < n; i += 4)
    {
        double x1 = x[i];
        double x2 = x[i+1];
        double x3 = x[i+2];
        double x4 = x[i+3];
        sum += 100 * pow(x1*x1 - x2, 2)
             + pow(x1 - 1, 2)
             + pow(x3 - 1, 2)
             + 90 * pow(x3*x3 - x4, 2)
             + 10.1 * (pow(x2 - 1, 2) + pow(x4 - 1, 2))
             + 19.8 * (x2 - 1) * (x4 - 1);
    }
    return sum;
}

double calc_powell_singular_function(const std::vector<double>& x, uint32_t n)
{
    double sum = 0.0;
    for (uint32_t i = 0; i + 3 < n; i += 4)
    {
        double x1 = x[i];
        double x2 = x[i+1];
        double x3 = x[i+2];
        double x4 = x[i+3];
        sum += pow(x1 + 10*x2, 2)
             + 5 * pow(x3 - x4, 2)
             + pow(x2 - 2*x3, 4)
             + 10 * pow(x1 - x4, 4);
    }
    return sum;
}

// ===== PUNKTY STARTOWE =====

std::vector<double> make_quadratic_x0(uint32_t n)
{
    return std::vector<double>(n, 5.0);
}

std::vector<double> make_woods_x0(uint32_t n)
{
    return std::vector<double>(n, -3.0);
}

std::vector<double> make_powell_x0(uint32_t n)
{
    return std::vector<double>(n, 3.0);
}
