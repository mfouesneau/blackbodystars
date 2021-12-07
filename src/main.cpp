/**
 * @file main.cpp
 * @brief Main example file
 * @version 0.1
 * @date 2021-11-23
 *
 */
#include <cmath>
#include <exception>
#include <fstream>
#include <iostream>
#include <ratio>
#include <sstream>
#include <stdexcept>
#include <string>
#include <xtensor/xarray.hpp>
#include <xtensor/xbuilder.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xview.hpp>
#include "rquantities.hpp"
#include "rapidcsv.hpp"
#include "blackbody.hpp"
#include "helpers.hpp"
#include "prettyprint.hpp"
#include <cpr/cpr.h>
#include "xml2json.hpp"
#include "rapidjson/document.h"
#include "votable.hpp"

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


std::string download_svo_filter(std::string id){
    cpr::Response r = cpr::Get(cpr::Url{"http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=Generic/Bessell_JHKLM.J"});
                               // cpr::Parameters{{"format", "ascii"}, {"id", "Generic/Bessell_JHKLM.J"}});
    std::cout << "\n URL: " << r.url << std::endl;
    std::cout << "\n Content: \n" << r.text << std::endl;
    return r.text;
}

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

std::vector<float> col = doc.GetColumn<float>("GALEX_FUV");
std::cout << col << "\n";
//std::cout << doc.GetColumnNames() << "\n";

// Seg fault in the following
// auto data = download_svo_filter("Generic/Bessell_JHKLM.J");

/*
rapidcsv::Document filt("/workspace/blackbodystars/data/passbands/Bessell_JHKLM.J",
                        rapidcsv::LabelParams(-1, -1),
                        rapidcsv::SeparatorParams(' ')
                        );
std::cout << filt.GetColumn<std::string>(0) << "\n";
std::vector<float> wavelength = doc.GetColumn<float>(0, rapidcsv::Converter<float>(rapidcsv::ConverterParams()));
//std::vector<double> transmission = doc.GetColumn<double>(1);

std::string xml_str, json_str;
std::ostringstream oss;
std::ifstream infile;
infile.open("data/passbands/GAIA.GAIA3.G.xml");
oss.str("");
oss << infile.rdbuf();
xml_str = oss.str().data();
infile.close();
json_str = xml2json(xml_str.c_str());

std::cout << json_str << "\n";

rapidjson::Document document;
document.Parse(json_str.c_str());
*/

votable::VOTable vot("data/passbands/GAIA.GAIA3.G.xml");
std::cout << vot.params["filterID"] 
            << " "
            << vot.get<double>("Wavelength")
            << "\n"
            << vot.get<double>(1)
            << "\n";

std::cout << "done.\n";
return 0;

}
