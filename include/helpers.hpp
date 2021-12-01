/**
 * @file helpers.hpp
 * @brief Some random tools to ease the project
 * @version 0.1
 * @date 2021-12-01
 *
 */
#pragma once
#include <iostream>
#include <vector>


/**
 * Display an `std::vector` object
 */
template <typename T>
std::ostream & operator<<(std::ostream &os,
                          const std::vector<T> &v)
{
    os << "[ ";
    for (const T &p : v){ os << p << " ";}
    os << "]\n";
    return os;
}


/**
 * Check if a vector contains an item
 */
bool contains(const std::vector<std::string>& v, const std::string& other){
    for (const auto& entry: v){
        if (entry.compare(other) == 0){
            return true;
        }
    }
    return false;
}