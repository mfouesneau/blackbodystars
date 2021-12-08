/**
 * @file io.hpp
 * @brief Tools to read and write filter libraries
 * @version 0.1
 *
 */
#pragma once
#include <cphot/filter.hpp>
#include <cphot/votable.hpp>
#include <cphot/rquantities.hpp>


namespace cphot {

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