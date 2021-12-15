/**
 * @defgroup LICKS Lick indices
 * @brief Lick indices calculations
 *
 * This package provides function to compute spectral indices.
 *
 * We provide a collection of many common indices. The Lick system of spectral
 * line indices is one of the most commonly used methods of determining ages and
 * metallicities of unresolved (integrated light) stellar populations.
 *
 * The calibration of the Lick/ IDS system is complicated because the original
 * Lick spectra were not flux calibrated, so there are usually systematic
 * effects due to differences in continuum shape.  Proper calibration involves
 * observing many of the original Lick/IDS standard stars and deriving offsets
 * to the standard system.
 *
 * references:
 * - Worthey G., Faber S. M., Gonzalez J. J., Burstein D., 1994, ApJS, 94, 687
 * - Worthey G., Ottaviani D. L., 1997, ApJS, 111, 377
 * - Puzia et al. 2002
 * - Zhang, Li & Han 2005, http://arxiv.org/abs/astro-ph/0508634v1
 *
 * @note
 * In Vazdekis et al. (2010), we propose a new Line Index System, hereafter LIS,
 * with three new spectral resolutions at which to measure the Lick indices.
 * Note that this new system should not be restricted to the Lick set of indices
 * in a flux calibrated system. In fact, LIS can be used for any index in the
 * literature (e.g., for the Rose (1984) indices), including newly defined
 * indices (e.g., Cervantes & Vazdekis 2009).
 * The LIS system is defined for 3 different spectral resolutions which are best
 * suited for the following astrophysical cases:
 * - LIS-5.0AA: globular clusters
 * - LIS-8.4AA: low and intermediate-mass galaxies
 * - LIS-14.0AA: massive galaxies
 * Conversions to transform the data from the Lick/IDS system to LIS can be found
 * with discussion of indices and information in Johansson, Thomas & Maraston (2010)
 * http://wwwmpa.mpa-garching.mpg.de/~jonasj/milesff/milesff.pdf
 *
 */
#include <cmath>
#include <cphot/rquantities.hpp>
#include <cphot/hardcoded_data/licks_data.hpp>

namespace cphot{

using DMatrix = xt::xarray<double, xt::layout_type::row_major>;

/**
 * @ingroup LICKS
 * @brief Adapt the resolution of the spectra to match the lick definitions.
 *
 * Lick definitions have different resolution elements as function of wavelength.
 * These definition are hard-coded in this function
 *
 * @param w
 *         wavelength definition in Angstrom
 * @param flux
 *         spectra to convert
 * @param fwhm0
 *         initial broadening in the spectra `fi`
 * @param sigma_floor
 *         minimal dispersion to consider
 * @return flux_red reduced spectra
 *
 */
DMatrix reduce_resolution(const DMatrix& w, const DMatrix& flux, double fwhm0, double sigma_floor){
    // all in AA
    const DMatrix w_lick_res = (4000., 4400., 4900., 5400., 6000.);   ///< Lick resolution anchor points in AA
    const DMatrix lick_res   = (11.5, 9.2, 8.4, 8.4, 9.8);              ///< FWHM in AA

    // Linear interpolation of lick_res over w
    // TODO: need to add extrapolation
    DMatrix res = xt::interp(w, w_lick_res, lick_res);

    // Compute width from fwhm
    double constant = 2. * std::sqrt(2. * std::log(2));     ///< constant that does the conversion fwhm --> sigma
    DMatrix lick_sigma = xt::sqrt((res * res - fwhm0 * fwhm0)) / constant;

    // Convolution by g=1/sqrt(2*pi*sigma^2) * exp(-r^2/(2*sigma^2))
    DMatrix flux_red = xt::zeros<double>(flux.shape);

    for (size_t i=0; i < lick_sigma.size(); ++i){
        double sigma = lick_sigma(i);
        double maxsigma = 3. * sigma;
        // sampling floor: min (0.2, sigma * 0.1)
        double delta = std::min(sigma_floor, sigma * 0.1);
        DMatrix delta_wj = xt::arange(-maxsigma, + maxsigma, delta);
        auto wj = delta_wj + w[i];
        auto fluxj = xt::interp(wj, w, flux, 0., 0.)
        flux_red[i] = xt::sum(fluxj * delta * xt::exp(-0.5 * xt::pow(delta_wj / sigma, 2)));
    }
    flux_red /= lick_sigma * const;
    return flux_red;
}

/**
 * @ingroup LICKS
 * @brief Define a Lick Index similarily to a Filter object
 */
class LickIndex{

    private:
        std::string name;             ///< name of the index
        double index_band_min;        ///< minimal wavelength of index interval
        double index_band_max;        ///< maximal wavelength of index interval
        double blue_continuum_min;    ///< minimal wavelength of blue continuum
        double blue_continuum_max;    ///< maximal wavelength of blue continuum
        double red_continuum_min;     ///< minimal wavelength of red continuum
        double red_continuum_max;     ///< maximal wavelength of red continuum
        QLength wavelength_unit;      ///< wavelength unit (usually Angstrom or nm)
        std::string description;      ///< description/notes

    public:
        LickIndex(const std::string& name, double index_band_min, double index_band_max,
                  double blue_continuum_min, double blue_continuum_max,
                  double red_continuum_min, double red_continuum_max,
                  const QLength & wavelength_unit, const std::string& description);
        LickIndex(const lickdata& data);

        std::string get_name() const;

        void info() const;
};

/**
 * @brief Construct a new Lick Index:: Lick Index object
 *
 * @param name                name of the index
 * @param index_band_min      minimal wavelength of index interval
 * @param index_band_max      maximal wavelength of index interval
 * @param blue_continuum_min  minimal wavelength of blue continuum
 * @param blue_continuum_max  maximal wavelength of blue continuum
 * @param red_continuum_min   minimal wavelength of red continuum
 * @param red_continuum_max   maximal wavelength of red continuum
 * @param wavelength_unit     wavelength unit (usually Angstrom or nm)
 * @param description         description/notes
 */
LickIndex::LickIndex(const std::string& name, double index_band_min, double index_band_max,
                     double blue_continuum_min, double blue_continuum_max,
                     double red_continuum_min, double red_continuum_max,
                     const QLength& wavelength_unit, const std::string& description)
            : name(name), index_band_min(index_band_min), index_band_max(index_band_max),
              blue_continuum_min(blue_continuum_min), blue_continuum_max(blue_continuum_max),
              red_continuum_min(red_continuum_min), red_continuum_max(red_continuum_max),
              wavelength_unit(wavelength_unit), description(description){}

/**
 * @brief Construct a new Lick Index:: Lick Index object
 *
 * @param data  data from the hardcoded data indices
 */
LickIndex::LickIndex(const lickdata& data)
            : name(data.name), index_band_min(data.index_band_min), index_band_max(data.index_band_max),
              blue_continuum_min(data.blue_continuum_min), blue_continuum_max(data.blue_continuum_max),
              red_continuum_min(data.red_continuum_min), red_continuum_max(data.red_continuum_max),
              wavelength_unit(data.wavelength_unit), description(data.description){}

/**
 * @brief Get the name object
 *
 * @return std::string
 */
std::string LickIndex::get_name() const {
    return this->name;
}


/**
 * @brief display some information about the index
 *
 */
void LickIndex::info(){
    double wave_conv = this->wavelength_unit.to(nm).value();

    std::cout << "Lick Index: " << this->name << "\n"
              << "  Index band: [" << this->index_band_min * wave_conv << ", " << this->index_band_max * wave_conv << "] nm" << "\n"
              << "  Blue continuum: [" << this->blue_continuum_min * wave_conv << ", " << this->blue_continuum_max * wave_conv << "] nm" << "\n"
              << "  Red continuum: [" << this->red_continuum_min * wave_conv << ", " << this->red_continuum_max * wave_conv << "] nm" << "\n";
}

/**
 * @brief Nice representation of Filter objects
 *
 * @param os   stream to output the representation
 * @param F    Filter object
 * @return std::ostream&  same as os
 */
std::ostream & operator<<(std::ostream &os,
                          LickIndex &index){
    os << "LickIndex: " << index.get_name()
       << "\n";
    return os;
}

/**
 * @ingroup LICKS
 * @brief Collection of Lick indices
 */
class LickLibrary{

};

}// namespace cphot