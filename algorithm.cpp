#include "algorithm.h"
#include <cmath>
#include <random>
#include <thread>

// ===== SEKWENCYJNE =====

std::pair<std::vector<double>, double>
perform_sequential_algorithm(const calc_function_t& calc_value,
                             std::vector<double> x0,
                             uint32_t n,
                             int a, int b)
{
    const uint32_t L = 30;
    double T = 500.0;
    const double alpha = 0.3;
    const double epsT = 0.1;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> U01(0.0, 1.0);

    double f_x0 = calc_value(x0, n);
    std::vector<double> xopt = x0;
    double f_opt = f_x0;

    while (T > epsT)
    {
        for (uint32_t k = 0; k < L; ++k)
        {
            std::vector<double> x_star(n);
            for (uint32_t i = 0; i < n; ++i)
                x_star[i] = a + U01(gen) * (b - a);

            double f_star = calc_value(x_star, n);

            if (f_star < f_x0 ||
                U01(gen) < std::exp((f_x0 - f_star) / T))
            {
                x0 = x_star;
                f_x0 = f_star;
                if (f_star < f_opt)
                {
                    xopt = x_star;
                    f_opt = f_star;
                }
            }
        }
        T *= (1.0 - alpha);
    }
    return {xopt, f_opt};
}

// ===== WIELOWĄTKOWE =====

std::pair<std::vector<double>, double>
perform_threaded_algorithm(const calc_function_t& calc_value,
                           uint32_t n,
                           int a, int b,
                           int num_threads)
{
    std::vector<std::thread> threads;
    std::vector<std::vector<double>> best_x(num_threads);
    std::vector<double> best_f(num_threads);

    for (int t = 0; t < num_threads; ++t)
    {
        threads.emplace_back([&, t]()
        {
            std::mt19937 gen(1000 + t);
            std::uniform_real_distribution<double> U(a, b);

            std::vector<double> x0(n);
            for (uint32_t i = 0; i < n; ++i)
                x0[i] = U(gen);

            auto result = perform_sequential_algorithm(
                calc_value, x0, n, a, b);

            best_x[t] = result.first;
            best_f[t] = result.second;
        });
    }

    for (auto& th : threads)
        th.join();

    int best = 0;
    for (int i = 1; i < num_threads; ++i)
        if (best_f[i] < best_f[best])
            best = i;

    return {best_x[best], best_f[best]};
}
