/**
 * @file cpr_test.cpp
 * @brief  testing CPR to download data from SVO
 * @version 0.1
 *
 */
#include <iostream>
#include "cphot/io.hpp"


int main(){
    cphot::Filter filt = cphot::download_svo_filter("2MASS/2MASS.H");
    filt.info();
}