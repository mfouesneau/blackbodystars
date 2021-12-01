/**
 * @file main.cpp
 * @brief Main example file
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
#include "rapidcsv.hpp"

/**
 * Blackbody as a flux distribution as function
 *  of wavelength, temperature and amplitude.
 *
 * @param lam:   wavelength in nm
 * @param amp:   dimensionless normalization factor
 * @param teff: temperature in Kelvins
 * @return evaluation of the blackbody radiation in flam units (erg/s/cm2/AA)
 *
 * Note that amp is alternatively represented as the angular size θ = R/d = sqrt(amp/pi)
 */
double bb_flux_function(double lam_nm,
                        double amp,
                        double teff_K){
    // Natural constants.
    double kB = 1.380649e-23;   // Unit("J/K")
    double c = 299792458.0;     // Unit("m/s")
    double h = 6.62607015e-34;  // Unit('m**2 * kg / s')
    double v = (amp * 2 * h * std::pow(c, 2) / (std::pow(lam_nm, 5) *
            (std::exp(h * c / (lam_nm * 1e-9 * kB * teff_K)) - 1)));
    std::cout << lam_nm << " " << amp << " " << teff_K << "\n";
    return v * 1e+38;  // flam = erg/s/cm2/AA
}

/**
 * Blackbody as a flux distribution as function
 *  of wavelength, temperature and amplitude.
 *
 * @param lam:   wavelength in nm
 * @param amp:   dimensionless normalization factor
 * @param teff: temperature in Kelvins
 * @return evaluation of the blackbody radiation in flam units (erg/s/cm2/AA)
 *
 * Note that amp is alternatively represented as the angular size θ = R/d = sqrt(amp/pi)
 *
 * References:
 *
 *  - Rybicki, G. B.; Lightman, A. P. (1979).
 *    _Radiative Processes in Astrophysics._ John Wiley & Sons.
 *    ISBN 0-471-82759-2.
 *
 * Planck's law of black-body radiation states that
 * \f{eqnarray*}{
 *      f_{\lambda }(T) &=& a\,\cdot\,{\frac {2h\nu ^{3}}{c^{2}}}{\frac {1}{e^{h\nu /kT}-1}}\\
 *                      &=& a\,\cdot\,\frac{2hc^2}{\lambda^5}{\frac {1}{e^{hc /\lambda kT}-1}},
 * \f}
 * where
 * - $f_\lambda(T)$ is the spectral flux density at thermal equilibrium of temperature $T$.
 * - $h$ is the Planck constant ($6.6260693\cdot10^{-34}$  J/s);
 * - $c$ is the speed of light in a vacuum ($299792458\cdot 10^6\mu m/s$);
 * - $k$ is the Boltzmann constant ($1.3806505\cdot 10^{-23}$ J/K);
 * - $\nu$ is the frequency of the electromagnetic radiation;
 * - $\lambda$ is the wavelength of the electromagnetic radiation;
 * - $T$ is the absolute temperature of the body;
 * - and $a$ is a dimensionless normalization factor.
 *
 * Note that $a$ can alternatively links to the angular size $\theta$
 * \f[
 *      \theta = \sqrt{\frac{a}{\pi}} = \frac{R}{d}
 * \f]
 * where $R$ is the radius, and $d$ is the distance to the star.
 *
 */
QSpectralFluxDensity bb_flux_function(QLength lam,
                                      Number amp,
                                      QTemperature teff){
    // Natural constants.
    double kB = (1.380649e-23 * joule / kelvin).getValue();
    double c = speed_of_light.Convert(meter / second);
    double h = (6.62607015e-34 * metre2 * kg / second).getValue();
    double lam_nm = lam.Convert(nanometre);
    double teff_K = teff.Convert(kelvin);
    double amp_ = amp.getValue();

    double v = (amp_ * 2 * h * std::pow(c, 2) / (std::pow(lam_nm, 5) *
            (std::exp(h * c / (lam_nm * 1e-9 * kB * teff_K)) - 1)));

    return v * 1e+38 * flam;
}


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

/**
 * Display an `std::map` object with <<
 */
std::ostream & operator <<(std::ostream &os,
                           const std::vector<std::string> &v)
{
    os << "[ ";
    for (const std::string &p : v){ os << p << " ";}
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


int main() {
    int A[3] {3, 5, 7};
    //auto a0, a1, a2] = A;
    for (auto val : A){
        std::cout << val << " ";
    }
    std::cout << "\n";

    example::example1();

    xt::xarray<double> arr1
  {{1.0, 2.0, 3.0},
   {2.0, 5.0, 7.0},
   {2.0, 5.0, 7.0}};

xt::xarray<double> arr2
  {5.0, 6.0, 7.0};

xt::xarray<double> res = xt::view(arr1, 1) + arr2;

std::cout << res << std::endl;

std::cout << bb_flux_function(500, 1., 5000) << std::endl;
std::cout << bb_flux_function(500e-9 * metre, 1., 5000 * kelvin).Convert(flam) << std::endl;

std::cout << (45_km).to(parsec) << "\n";

std::cout << (1 * lsun) << "\n";

std::cout << "45 km = " << (45_km).to(parsec) << " pc \n";

std::cout << "1 Lsun " << (1 * lsun) << " W \n";  // default definition

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
