/**
 * @file main.cpp
 * @brief Main example file
 * @version 0.1
 * @date 2021-11-23
 *
 */
#include <iostream>
#include <exception>
#include <cmath>
#include <ratio>
#include <stdexcept>
#include <xtensor/xarray.hpp>
#include <xtensor/xbuilder.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xview.hpp>
#include "rquantities.hpp"
#include "rapidcsv.hpp"
#include "blackbody.hpp"
#include "helpers.hpp"
#include "prettyprint.hpp"

using DMatrix = xt::xarray<double, xt::layout_type::row_major>;

class Filter {

    private:
        static constexpr double c = speed_of_light.to(angstrom / second);
        static constexpr double h = 6.62607015e-27;  // erg/s

        DMatrix wavelength_nm;
        DMatrix transmission;
        std::string name = "";
        std::string dtype = "photon";

    public:
        Filter(const DMatrix& wavelength,
               const DMatrix& transmission,
               const QLength& wavelength_unit,
               const std::string dtype,
               const std::string name){
                    double convfac = wavelength_unit.to(nanometre);
                    this->wavelength_nm = convfac * wavelength;
                    this->name = name;
                    this->transmission = transmission;

                    if ((dtype.compare("photon") == 0) ||
                        (dtype.compare("energy"))){
                        this->dtype = dtype;
                    } else {
                        throw std::runtime_error("only photon and energy allowed");
                    }
                }

};

int main() {

rapidcsv::Document doc("data/blackbody-stars-clean.csv");

std::vector<std::string> ignore = {
    "SDSSName", "R.A.(J2000)", "Decl.(J2000)", "Teff",
    "Teff_error", "amp", "amp_error", "theta", "theta_error",
    "chi^2/dof", "source_id",
    };

std::vector<std::string> columns;
for (const auto& p: doc.GetColumnNames()){
    if (! contains(ignore, p)){ columns.push_back(p); }
}
std::cout << columns << "\n";

// std::vector<float> col = doc.GetColumn<float>("GALEX_FUV");
// std::cout << doc.GetColumnNames() << "\n";
std::cout << "done.\n";

}
