/**
 * @brief: Reading VOTABLEs
 *
 * version: 0.1a
 */
#pragma once
#include "xml2json.hpp"
#include "rapidjson/document.h"
#include <string>
#include <map>


template <typename T>
T parseString(const char * txt){
    T d  = std::atof(txt);
    return d;
}

template <typename T>
T parseString(const std::string& txt){
    T d  = std::atof(txt.c_str());
    return d;
}

namespace votable {

    class VOTable {

        public:
            VOTable(const std::string & input_filename);
            std::string version;
            std::map<std::string, std::string> param;
            std::map<std::string, std::string> param_type;



        private:
            rapidjson::Document document;  /** Storing the JSON document */
            std::string get_version(){return document["VOTABLE"]["@version"].GetString();}
            void testing();

    };


/**
 * Constructor
 *
 * @param input_filename: XML file to parse
 */
VOTable::VOTable(const std::string & input_filename){
    std::ostringstream oss;
    std::ifstream infile;
    infile.open(input_filename);
    oss.str("");
    oss << infile.rdbuf();
    infile.close();
    std::string json_str = xml2json(oss.str().data());
    this->document.Parse(json_str.c_str());
    this->version = this->get_version();
    this->testing();
}

void VOTable::testing(){
    const rapidjson::Value& where = this->document["VOTABLE"]["RESOURCE"]["TABLE"]["PARAM"];
    for (auto& v : where.GetArray()){
        if (std::string(v["@datatype"].GetString()).compare("float") == 0){
            //double d  = std::atof(v["@value"].GetString());
            std::cout << v["@name"].GetString() << ": " 
                    << parseString<double>(v["@value"].GetString())
                    // << " (" << v["@datatype"].GetString() << ")"
                    << "\n";
        } else {
            std::cout << v["@name"].GetString() << ": " 
                    << v["@value"].GetString() 
                    << " (" << v["@datatype"].GetString() << ")"
                    << "\n";
        }
    }
}


} // namespace votable