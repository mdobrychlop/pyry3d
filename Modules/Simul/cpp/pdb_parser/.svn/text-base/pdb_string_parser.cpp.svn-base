
#include <string>
#include <algorithm>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include "../pyry3d_cpp.hpp"
#include "pdb_input_parser.hpp"

// trims string from left (whitespaces)
std::string& ltrim(std::string& s) {
	s.erase(s.begin(), find_if(s.begin(), s.end(), not1(std::ptr_fun<int, int>(isspace))));
	return s;
}

// trims from right (whitespaces)
std::string& rtrim(std::string& s) {
	s.erase(find_if(s.rbegin(), s.rend(), not1(std::ptr_fun<int, int>(isspace))).base(), s.end());
	return s;
}

// // trims from both sides (whitespaces)
std::string& trim(std::string& s) {
	return ltrim(rtrim(s));
}

std::vector<pdb_line> parse_first_file(const std::string& pdb_string) {
	/* interprets pdb_string to cpp structs*/
	std::vector<pdb_line> already_read;

	if (!pdb_string.size()) 
		std::cerr<< "No pdb_string!!\n";

	std::stringstream ifs(pdb_string);

	std::string line, record_name, atom_ids, residue_name, atom_name;
	std::string chain_identifiers, residue_numbers;
	std::string coordsx, coordsy, coordsz;
	std::string occupancy, temperature_factor, elem_symbol;
	
	while (getline(ifs, line)) {
		
		// HAHA, MAGICZNE 78, LAURA, PROSZĘ ZMIEŃ TO
		if (line.size() < 78) continue;

		record_name = line.substr(0, 6);
		atom_ids = line.substr(6, 5);
		atom_name = line.substr(12, 4);
		residue_name = line.substr(17, 3);

		chain_identifiers = line.substr(21, 1);
		residue_numbers = line.substr(22, 4); // 23 - 26

		coordsx = line.substr(30, 8); // 31 - 38
		coordsy = line.substr(38, 8);
		coordsz = line.substr(46, 8); // 47 - 54

		occupancy = line.substr(54, 6); // 55 - 60
		temperature_factor = line.substr(60, 6); // 61 - 66 Default = 0.0
		elem_symbol = line.substr(76, 2); // 77 - 78 right-justified

		trim(record_name);
		trim(atom_ids);
		trim(residue_name);

		trim(atom_name);
		trim(chain_identifiers);
		trim(residue_numbers);

		trim(coordsx);
		trim(coordsy);
		trim(coordsz);

		trim(occupancy);
		trim(temperature_factor);
		trim(elem_symbol);

		pdb_line line_obj(record_name, atom_ids, atom_name, residue_name,
			chain_identifiers, residue_numbers,
			coordsx, coordsy, coordsz,
			occupancy, temperature_factor, elem_symbol);

        already_read.push_back(line_obj);
	}
	return already_read;
}

