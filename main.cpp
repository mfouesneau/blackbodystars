#include <iostream>
#include "example.hpp"
#include <iostream>
#include <cmath>
#include "xtensor/xarray.hpp"
#include "xtensor/xio.hpp"
#include "xtensor/xview.hpp"

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
    return v * 1e+38;  // flam = erg/s/cm2/AA  
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
    
}
