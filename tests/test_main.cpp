/**
 * @file test_main.cpp
 * @brief Unittest on Main example file
 * @version 0.1
 * @date 2021-11-23
 *
 */
#include <iostream>
#include "example.hpp"
#include <exception>
#include <iostream>
#include <cmath>
#include <ratio>
#include <stdexcept>
#include <xtensor/xarray.hpp>
#include <xtensor/xbuilder.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xview.hpp>
#include "rquantities.hpp"

void test_array(){
    int A[3] {3, 5, 7};
    for (auto val : A){
        std::cout << val << " ";
    }
    std::cout << "\n";
}

void test_example1(){
    example::example1();
}

void test_xtensor(){
    xt::xarray<double> arr1
    {{1.0, 2.0, 3.0},
    {2.0, 5.0, 7.0},
    {2.0, 5.0, 7.0}};

    xt::xarray<double> arr2
    {5.0, 6.0, 7.0};

    xt::xarray<double> res = xt::view(arr1, 1) + arr2;

    std::cout << res << std::endl;
}

void test_units(){
    std::cout << (45_km).to(parsec) << "\n";
    std::cout << (1 * lsun) << "\n";
}

int main() {
    test_array();
    test_example1();
    test_xtensor();
    test_units();
    // std::cout << bb_flux_function(500, 1., 5000) << std::endl;
    // std::cout << bb_flux_function(500e-9 * metre, 1., 5000 * kelvin).Convert(flam) << std::endl;
    return 0;
}
