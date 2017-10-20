#ifndef PDB_STRING_PARSER_PYRY3D
#define PDB_STRING_PARSER_PYRY3D

#include <string>
#include "pdb_input_parser.hpp"

std::string& ltrim(std::string& s);
std::string& rtrim(std::string& s);
std::string& trim(std::string& s);

std::vector<pdb_line> parse_first_file(const std::string& pdb_string);

void parse_file(const std::string& pdb_string, std::vector<point>& coords);

std::string coord_to_str(float value);

#endif
