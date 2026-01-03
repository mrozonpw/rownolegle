#include <iostream>
#include <chrono>
#include "algorithm.h"
#include "functions.h"

int main()
{
    uint32_t n = 200000;

    // ===== SEKWENCYJNIE =====
    auto x0 = make_quadratic_x0(n);
    auto t1 = std::chrono::high_resolution_clock::now();
    auto seq = perform_sequential_algorithm(
        calc_quadratic_function, x0, n, -5, 5);
    auto t2 = std::chrono::high_resolution_clock::now();

    auto T1 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cout << "SEQ time: " << T1 << " ms, f = " << seq.second << "\n";

    // ===== WĄTKI =====
    for (int p : {2, 4, 8})
    {
        auto t3 = std::chrono::high_resolution_clock::now();
        auto par = perform_threaded_algorithm(
            calc_quadratic_function, n, -5, 5, p);
        auto t4 = std::chrono::high_resolution_clock::now();

        auto Tp = std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count();
        std::cout << p << " threads: " << Tp
                  << " ms, speedup = "
                  << double(T1) / Tp
                  << "\n";
    }

    return 0;
}
