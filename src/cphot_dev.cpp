/**
 * @file cphot_dev.cpp
 * @brief  developping cphot
 * @version 0.1
 *
 */
#include <iostream>
#include <cphot/io.hpp>
#include <cphot/votable.hpp>


/**
 * @brief Interface to Vega reference data
 *
 *   Class that handles vega spectrum and references.  This class know where to
 *   find the Vega synthetic spectrum (Bohlin 2007) in order to compute fluxes
 *   and magnitudes in given filters
 *
 *   Attributes
 *   ----------
 *   source: str
 *       filename of the vega library
 *   data: SimpleTable
 *       data table
 *   units: tuple
 *       detected units from file header
 *   wavelength: array
 *       wavelength (with units when found)
 *   flux: array
 *       flux(wavelength) values (with units when provided)
 *
 **/
class Vega {
    public:
        Vega(const std::string& source);
        Vega();

};

Vega::Vega() {
}

Vega::Vega(const std::string& source) {
}


int main(){
    cphot::Filter filt = cphot::download_svo_filter("2MASS/2MASS.H");
    filt.info();

    votable::VOTable vega("vega.vot");
    std::cout << vega << std::endl;

}