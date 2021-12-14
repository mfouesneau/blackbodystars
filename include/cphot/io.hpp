/**
 * @defgroup IO Input/Output
 * @brief Tools to read and write filter libraries
 * @version 0.1
 *
 */
#pragma once
#include "filter.hpp"
#include "votable.hpp"
#include "rquantities.hpp"
#include <cpr/cpr.h>


namespace cphot {

/**
 * @ingroup IO
 * @brief Get the filter object from VOTable object
 *
 * @param vot  votable document
 * @return Filter object
 */
Filter get_filter(votable::VOTable& vot){
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

/**
 * @ingroup IO
 * @brief Get the filter object from a VOTable file
 *
 * @param vot_filename  path to the xml file
 * @return Filter object
 */
Filter get_filter(const std::string& vot_filename){
    // Read VOTable
    votable::VOTable vot(vot_filename);
    return get_filter(vot);
}

/**
 * @ingroup IO
 * @brief main interface to SVO data requests
 *
 * Query the <a link=http://svo2.cab.inta-csic.es/theory/fps/>SVO filter profile
 * service</a> and return the filter object
 * (http://svo2.cab.inta-csic.es/theory/fps)
 *
 * @param id              passband id
 * @return std::string    data content
 */
Filter download_svo_filter(const std::string & id){
    cpr::Response r = cpr::Get(cpr::Url{"http://svo2.cab.inta-csic.es/theory/fps/fps.php"},
                               cpr::Parameters{{"ID", id}});
    votable::VOTable vot;
    vot.from_content(r.text);
    return cphot::get_filter(vot);
}

}; // namespace cphot