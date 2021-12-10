/**
 * @file cphot_dev.cpp
 * @brief  developping cphot
 * @version 0.1
 *
 */
#include <iostream>
#include <cphot/io.hpp>
#include <cphot/rquantities.hpp>
#include <cphot/votable.hpp>
#include <cphot/vega.hpp>



int main(){
    std::string filter_id = "2MASS/2MASS.H";
    cphot::Filter filt = cphot::download_svo_filter(filter_id);
    filt.info();

    cphot::Vega v2 = cphot::Vega(
        cphot_vega::wavelength_nm,
        cphot_vega::flux_flam,
        nm, flam);

    double flux_flam_v2 = filt.get_flux(v2.get_wavelength(nm), v2.get_flux(flam), nm, flam).to(flam);
    std::cout << "Vega zero points for filter: " << filter_id << "\n"
              <<  flux_flam_v2 << " flam\n"
              << -2.5 * std::log10(flux_flam_v2) << " mag\n";

}