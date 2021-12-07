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
#include "prettyprint.hpp"


/**
 * @brief Parse a char* into a numeric value.
 * 
 * @tparam T   the output type 
 * @param txt  the string to parse
 * @return T   the value as a T type
 */
template <typename T>
T parseString(const char * txt){
    T d  = std::atof(txt);
    return d;
}

/**
 * @brief Parse a std::string into a numeric value.
 * 
 * @tparam T   the output type 
 * @param txt  the string to parse
 * @return T   the value as a T type
 */
template <typename T>
T parseString(const std::string& txt){
    T d  = std::atof(txt.c_str());
    return d;
}

/**
 * @brief Store Table Parameter Attributes
 * 
 */
struct Param {
    std::string name;
    std::string datatype;
    std::string value;
    std::string ucd;
    std::string utype;
    std::string unit;
    std::string description;
};

/**
 * @brief Store Table Field and data
 * 
 */
struct Field {
    std::string name;
    std::string datatype;
    std::string ucd;
    std::string utype;
    std::string unit;
    std::string description;
    std::vector<std::string> data;
};


/**
 * Display an `std::vector` object
 */
std::ostream & operator<<(std::ostream &os,
                          const Param &v)
{
    os << "Param: " << v.name << "= "
       << v.value  << " [" << v.unit << "]\n"
       << "     DTYPE=" << v.datatype << "\n"
       << "     UCD=" << v.ucd << "\n"
       << "     DESCRIPTION=" << v.description << "\n";
    return os;
}


/**
 * Display an `std::vector` object
 */
std::ostream & operator<<(std::ostream &os,
                          const Field &v)
{
    os << "Field: " << v.name << "length="
       << v.data.size()  << " [" << v.unit << "]\n"
       << "     DTYPE=" << v.datatype << "\n"
       << "     UCD=" << v.ucd << "\n"
       << "     DESCRIPTION=" << v.description << "\n";
    return os;
}



namespace votable {

    class VOTable {

        public:
            VOTable(const std::string & input_filename);
            std::string version;
            std::map<std::string, Param> params;  ///< Table parameters
            std::map<std::string, Field> fields;  ///< Table parameters

        private:
            rapidjson::Document document;  /** Storing the JSON document */
            std::string get_version(){return document["VOTABLE"]["@version"].GetString();}
            void testing();
            void setup_params();
            void setup_fields();

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
    this->setup_params();
    this->setup_fields();
    this->version = this->get_version();

    this->testing();
}

void VOTable::setup_params(){

    // read in the document PARAMs
    const rapidjson::Value& where = this->document["VOTABLE"]["RESOURCE"]["TABLE"]["PARAM"];
    for (const auto& v : where.GetArray()){
        Param p;
        auto name = v["@name"].GetString();
        p.name = name;
        p.datatype = v["@datatype"].GetString();
        p.value = v["@value"].GetString();
        if (v.HasMember("@ucd")) {p.ucd = v["@ucd"].GetString();}
        if (v.HasMember("@ytype")) {p.utype = v["@utype"].GetString();}
        if (v.HasMember("@unit")) {p.unit = v["@unit"].GetString();}
        if (v.HasMember("@DESCRIPTION")) {p.description = v["@DESCRIPTION"].GetString(); }
        this->params[name] = p;
    }
}

void VOTable::setup_fields(){ 

    // Read in FIELD information
    const rapidjson::Value& where = this->document["VOTABLE"]["RESOURCE"]["TABLE"]["FIELD"];
    for (const auto& v : where.GetArray()){
        Field p;
        auto name = v["@name"].GetString();
        p.name = name;
        p.datatype = v["@datatype"].GetString();
        if (v.HasMember("@ucd")) {p.ucd = v["@ucd"].GetString();}
        if (v.HasMember("@ytype")) {p.utype = v["@utype"].GetString();}
        if (v.HasMember("@unit")) {p.unit = v["@unit"].GetString();}
        if (v.HasMember("@DESCRIPTION")) {p.description = v["@DESCRIPTION"].GetString(); }
        this->fields[name] = p;
    }
    std::cout << this->fields << "\n";
}

void VOTable::testing(){
    std::cout << this->params["filterID"] 
              << " "
              << this->params["DetectorType"] 
              << " "
              << this->params["MagSys"] 
              << "\n";
}


} // namespace votable