/**
 * @file cpr_test.cpp
 * @brief  testing CPR to download data from SVO
 * @version 0.1
 * 
 */
#include <iostream>
#include <cpr/cpr.h>

/**
 * @brief main interface to SVO data requests
 * 
 * @param id              passband id
 * @return std::string    data content
 */
std::string download_svo_filter(std::string id){
    cpr::Response r = cpr::Get(cpr::Url{"http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=votable&id=Generic/Bessell_JHKLM.J"});
                               // cpr::Parameters{{"format", "ascii"}, {"id", "Generic/Bessell_JHKLM.J"}});
    std::cout << "\n URL: " << r.url << std::endl;
    std::cout << "\n Content: \n" << r.text << std::endl;
    return r.text;
}


int main(){

   cpr::Parameters{{"format", "ascii"}, {"id", "Generic/Bessell_JHKLM.J"}}; 

}