#include <cphot/filter.hpp>
#include <cphot/library.hpp>

int main(){

    std::string filter_name ("GaiaDR2_BP");
    std::string filename ("new_filters.hd5");
    auto tmp = cphot::get_filter_from_hdf5_library(filename, filter_name);

    tmp.info();

    auto lib = cphot::HDF5Library(filename);

    std::cout << lib << "\n";

    for (auto & c : lib.get_content()) {
       lib.load_filter(c).info();
    }
    std::cout << lib.find("gaia", false) << "\n";

    return 0;
}