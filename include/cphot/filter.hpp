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

        //! wavelength of the filter stored in nm
        DMatrix wavelength_nm;
        //! transmission of the passband
        DMatrix transmission;
        //! name of the filter
        std::string name = "";
        //! type, either photon or energy
        std::string dtype = "photon";
        //! units of the wavelength (nm by construction)
        QLength wavelength_unit = nm;
        //! Central wavelength in nm
        double cl;
        //! pivot wavelength in nm
        double lpivot;
        //! minimum wavelength in nm
        double lmin;
        //! maximum wavelength in nm
        double lmax;
        //! norm of the passband
        double norm;
        //! effective width in nm
        double width;
        //! full width at half maximum in nm
        double fwhm;
        //! Photon distribution based effective wavelength.
        double lphot = 0;
        //! Effective wavelength
        double leff = 0;
        //! Internal int λ * transmission * dλ
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
        QLength get_leff();
        QLength get_lphot();
        QLength get_fwhm();
        QLength get_width();
        QLength get_norm();
        QLength get_lmax();
        QLength get_lmin();
        QLength get_lpivot();
        QLength get_cl();
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
    this->wavelength_unit = nm;
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
 * @brief  Central wavelength
 *
 * \f[ \lambda_{cl} = \frac{\int \lambda T(\lambda) d\lambda}{\int T(\lambda) d\lambda}\f]
 *
 * @return central wavelength in nm
 */
QLength Filter::get_cl(){ return this->cl * this->wavelength_unit;}

/**
 * @brief  Pivot wavelength in nm
 *
 * if photon detector:
 * \f[ \lambda_p^2 = \frac{\int \lambda T(\lambda) d\lambda}{\int T(\lambda) d\lambda / \lambda} \f]
 *
 * if energy:
 * \f[ \lambda_p^2 = \frac{\int T(\lambda) d\lambda}{\int T(\lambda) d\lambda / \lambda^2} \f]
 *
 * @return pivot wavelength in nm
 */
QLength Filter::get_lpivot(){ return this->lpivot * this->wavelength_unit;}

/**
 * @brief the first λ value with a transmission at least 1% of maximum transmission
 *
 * @return min wavelength in nm
 */
QLength Filter::get_lmin(){ return this->lmin * this->wavelength_unit;}

/**
 * @brief the last λ value with a transmission at least 1% of maximum transmission
 *
 * @return max wavelength in nm
 */
QLength Filter::get_lmax(){ return this->lmax * this->wavelength_unit;}

/**
 * @brief the norm of the passband
 *
 * \f[ norm = \int T(\lambda) d\lambda \f]
 *
 * @return norm
 */
QLength Filter::get_norm(){ return this->norm; }

/**
 * @brief  Effective width
 *
 * \f[ width = \frac{\int T(\lambda) d\lambda}{\max(T(\lambda))} \f]
 *
 * @return width in nm
 */
QLength Filter::get_width(){ return this->width * this->wavelength_unit;}

/**
 * @brief the difference between the two wavelengths for which filter
 * transmission is half maximum.
 *
 * ..note::
 *      This calculation is not exact but rounded to the nearest passband
 *      data points
 *
 * @return fwhm in nm
 */
QLength Filter::get_fwhm(){ return this->fwhm * this->wavelength_unit;}

/**
 * @brief Photon distribution based effective wavelength.
 *
 * Defined as
 * \f[ \lambda_{phot} = \frac{\int\lambda^2 T(\lambda) Vega(\lambda) d\lambda }{\int\lambda T(\lambda) Vega(\lambda) d\lambda} \f]
 *
 * which we calculate as
 * \f[ \lambda_{phot} = \frac{get\_flux(\lambda Vega(\lambda))}{get\_flux(Vega(\lambda))} \f]
 *
 * @return QLength
 */
QLength Filter::get_lphot(){ return this->lphot * this->wavelength_unit;}

/**
 * @brief Effective wavelength
 *
 * \f[ \lambda_{eff} = \frac{\int \lambda T(\lambda) Vega(\lambda) d\lambda)}{\int T(\lambda) Vega(\lambda) d\lambda)} \f]
 *
 * @return Effective wavelenth
 */
QLength Filter::get_leff(){ return this->leff * this->wavelength_unit;}


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