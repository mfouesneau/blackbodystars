/**
 * @file filter.hpp
 * @brief A Filter class with unit awareness
 * @version 0.1
 */
 #pragma once
#include "rquantities.hpp"
#include <cmath>
#include <exception>
#include <iostream>
#include <regex>
#include <stdexcept>
#include <string>
#include <xtensor/xadapt.hpp>
#include <xtensor/xarray.hpp>

namespace cphot {

using DMatrix = xt::xarray<double, xt::layout_type::row_major>;

/**
 * @brief Unit Aware Filter.
 * input spectra and output values have units to avoid mis-interpretation.  
 * 
 * Note the usual (non SI) units of flux definitions:
 *       flam     = erg/s/cm**2/AA
 *       fnu      = erg/s/cm**2/Hz
 *       photflam = photon/s/cm**2/AA
 *       photnu   = photon/s/cm**2/Hz
 *
 * Define a filter by its name, wavelength and transmission The type of
 * detector (energy or photon counter) can be specified for adapting
 * calculations. (default: photon)
 */
class Filter {

    private:
        static constexpr double c = speed_of_light.to(angstrom / second);
        static constexpr double h = 6.62607015e-27;  // erg/s

        DMatrix wavelength_nm;
        DMatrix transmission;
        std::string name = "";
        std::string dtype = "photon";

        double cl;
        double lpivot;
        double lmin;
        double lmax;
        double norm;
        double width;
        double fwhm;
        double lphot = 0;
        double leff = 0;
        double lT;

        void calculate_sed_independent_properties();

    public:
        Filter(const DMatrix& wavelength,
               const DMatrix& transmission,
               const QLength& wavelength_unit,
               const std::string dtype,
               const std::string name); 
        std::string get_name(){ return this->name;}
        void info();
};

/**
 * @brief Construct a new Filter:: Filter object
 * 
 * @param wavelength       wavelength definition
 * @param transmission     transmission on the wavelength
 * @param wavelength_unit  units of the wavelength definition
 * @param dtype            "photon" or "energy"
 * @param name             name of the passband
 * @throw std::runtime_error if detector type is invalid
 */
Filter::Filter(const DMatrix& wavelength,
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
    this->calculate_sed_independent_properties();
}

/**
 * @brief Nice representation of Filter objects
 * 
 * @param os   stream to output the representation
 * @param F    Filter object
 * @return std::ostream&  same as os
 */
std::ostream & operator<<(std::ostream &os,
                          Filter &F){
    
    os << "Filter: " << F.get_name() 
       << "\n";
    return os;
}


/**
 * @brief Calculate the various standard properties of a given filter.
 * 
 * These properties are e.g., fwhm, pivot wavelength.
 * Those that do not require to consider an SED such as Vega.
 */
void Filter::calculate_sed_independent_properties(){

    // Calculate Filter properties
    const auto& wavelength_nm = this->wavelength_nm;
    const auto& transmission = this->transmission;

    size_t n_points = transmission.size();
    double transmission_max = xt::amax(transmission)[0];
    double transmission_max_100th = transmission_max / 100.;

    auto norm = xt::trapz(transmission, wavelength_nm)[0];
    auto _lT = xt::trapz(wavelength_nm * transmission, wavelength_nm)[0];
    auto _cl = norm > 0 ? _lT / norm : 0.;
    this->cl = _cl;
    this->norm = norm;
    this->lT = _lT;
    double lpivot2 = 0.;
    if (this->dtype.compare("photon") == 0){
        lpivot2 = _lT / trapz(transmission / wavelength_nm, wavelength_nm)[0];
    } else {
        lpivot2 = norm / trapz(transmission / xt::square(wavelength_nm), wavelength_nm)[0];
    }
    this->lpivot = std::sqrt(lpivot2);

    // the last value with a transmission at least 1% of maximum transmission
    double lmax = wavelength_nm[0];
    // the first value with a transmission at least 1% of maximum transmission
    double lmin = wavelength_nm[n_points - 1];
    for (size_t i=0; i < transmission.size(); ++i){
        if (transmission[i] > transmission_max_100th){
            lmax = std::max(lmax, wavelength_nm[i]);
            lmin = std::min(lmin, wavelength_nm[i]);
        }
    }
    this->lmin = lmin;
    this->lmax = lmax;

    // Effective width
    // Equivalent to the horizontal size of a rectangle with height equal
    // to maximum transmission and with the same area that the one covered by
    // the filter transmission curve.
    // W = int(T dlamb) / max(T)
    this->width = (norm / xt::amax(transmission)[0]);


    // FWHM
    // the difference between the two wavelengths for which filter transmission is
    // half maximum
    //
    // ..note:: 
    //      This calculation is not exact but rounded to the nearest passband data
    //      points
    double first = wavelength_nm[0];
    double last = wavelength_nm[-1];
    double thresh = transmission_max * 0.5;
    for (size_t i=0; i < wavelength_nm.size() - 1; ++i){
        if((transmission[i+1] > thresh) and (transmission[i] <= thresh)){ 
            first = wavelength_nm[i]; 
            break;
        }
    }
    for (size_t i=wavelength_nm.size(); i > 1; --i){
        if((transmission[i-1] > thresh) and (transmission[i] <= thresh)){ 
            last = wavelength_nm[i]; 
            break;
        }
    }
    this->fwhm = last - first;

    // leff = int (lamb * T * Vega dlamb) / int(T * Vega dlamb)
    // TODO: need vega 

    // lphot = int(lamb ** 2 * T * Vega dlamb) / int(lamb * T * Vega dlamb)
    // which we calculate as lphot = get_flux(lamb * vega) / get_flux(vega)
}

/**
 * @brief Display some information on cout
 */
void Filter::info(){
    size_t n_points = this->transmission.size();
    std::cout << "Filter Object information:\n"
            << "    name:                " << this->name << "\n"
            << "    detector type:       " << this->dtype << "\n"
            << "    wavelength units:    " << " nm  (internally set)" << "\n"
            << "    central wavelength:  " << this->cl  << " nm" << "\n"
            << "    pivot wavelength:    " << this->lpivot << " nm" << "\n"
            << "    minimum wavelength:  " << this->lmin << " nm" << "\n"
            << "    maximum wavelength:  " << this->lmax << " nm" << "\n"
            << "    norm:                " << this->norm << "\n"
            << "    effective width:     " << this->width << " nm" << "\n"
            << "    fullwidth half-max:  " << this->fwhm  << " nm" << "\n"
            << "    definition contains " << n_points << " points" << "\n";
}


/**
 * @brief Get the filter object from a VOTable file
 * 
 * @param vot_filename  path to the xml file
 * @return Filter object
 */
Filter get_filter(const std::string& vot_filename){
    // Read VOTable
    votable::VOTable vot(vot_filename);
    // Extract name
    std::string filter_name = std::regex_replace(
                                vot.params["filterID"].value, 
                                std::regex("/"),
                                "_"
                                );
    //extract detector type
    std::string detector_type = (parseString<int>(vot.params["DetectorType"].value) == 0) ? 
                                "energy" : "photon";

    // extract wavelength and units
    const std::vector<std::string> AA_str {"Angstrom", "AA", "angstrom"};
    const std::vector<std::string> nm_str {"Nanometer", "nanometer", "nm"};
    const auto & wave = vot.get<double>("Wavelength");
    QLength wavelength_unit;
    if (contains(AA_str, wave.unit)) { wavelength_unit = angstrom;}
    else if (contains(nm_str, wave.unit)) { wavelength_unit = nm;}

    // extract transmission
    const auto & transmit = vot.get<double>("Transmission");

    // convert to Filter inputs
    std::vector<std::size_t> shape = { wave.data.size() };
    DMatrix xt_wave = xt::adapt(wave.data, shape);
    DMatrix xt_transmit = xt::adapt(transmit.data, shape);

    return Filter(xt_wave, xt_transmit, 
                  wavelength_unit, detector_type, 
                  filter_name);
}


}; // namespace cphot