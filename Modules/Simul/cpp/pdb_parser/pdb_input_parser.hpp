#ifndef PDB_INPUT_PARSER_PYRY3D
#define PDB_INPUT_PARSER_PYRY3D

#include "../extern.hpp"
#include "../pyry3d_cpp.hpp"
#include <string>
#include <vector>

struct pdb_line {
	std::string record_name;
	int atom_ids;
	std::string residue_name;
	std::string atom_name;
	std::string chain_identifiers;
	int residue_numbers;
	float coordsx;
	float coordsy;
	float coordsz;
	float occupancy;
	float temperature_factor;
	std::string elem_symbol;

	pdb_line(std::string record_name, std::string atom_ids, std::string atom_name,
		std::string residue_name,
		std::string chain_identifiers, std::string residue_numbers,
		std::string coordsx, std::string coordsy, std::string coordsz,
		std::string occupancy, std::string temperature_factor, std::string elem_symbol);

};

int parse_pdb(const std::string& pdb_string);
void fill_all_array(const std::vector<pdb_line>& already_read);
void fill_res_to_ind_array();

#endif
