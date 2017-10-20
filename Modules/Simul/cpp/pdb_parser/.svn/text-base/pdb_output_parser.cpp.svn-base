#include <string>
#include <cstring>
#include <algorithm>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "pdb_input_parser.hpp"
#include "pdb_string_parser.hpp"
#include "../logger/logger.hpp"

int write_to_file(int num, const std::string& file_name, int atoms_number) {
	
	/* write complex from cpp module to pdb file */
	
	logger.write_to_logfile(num);

	char* cstr = new char[file_name.length() + 1];
	std::strcpy (cstr, file_name.c_str());
	
	FILE * file = fopen(cstr, "w");
	
	std::string prev_chain = "AA"; // there is no AA chain identifier
	const char* record_name = "ATOM";
	std::string s_atom_id, s_atom_name;
	std::string s_atom_coordx, s_atom_coordy, s_atom_coordz;
	std::string s_atom_occupancies, s_atom_temp_factors, s_atom_elem_symbols;
	std::string s_residue_name, s_residue_numbers, s_chain_identifier;


	char* at_name = new char[5];
	char* res_name = new char[6];
	char* chain_id = new char[2];
	char* at_symbols = new char[2];
	char* pre_chain_cstring = new char[3];
	
	for (int i = 0; i < atoms_number; i++) {
		
		std::strcpy (chain_id, chain_identifiers[i].c_str());
		std::strcpy(pre_chain_cstring, prev_chain.c_str());
		if (prev_chain != chain_id && prev_chain != "AA") {
			fprintf(file, "TER\n");
		}
		
		std::strcpy (at_name, atom_names[i].c_str());
		std::strcpy (res_name, residue_name[i].c_str());
		std::strcpy (at_symbols, atom_elem_symbols[i].c_str());
			
		prev_chain = chain_id;

		fprintf(file, "%-6s", record_name); // 1-6
		fprintf(file, "%5d", atom_ids[i]); // 7-11
		fprintf(file, "  ");
		fprintf(file, "%-4s", at_name); //13-16
		fprintf(file, "%3s", res_name); //18-20
		fprintf(file, "%2s", chain_id); //22
		fprintf(file, "%4d", residue_numbers[i]); //23-26
		
		fprintf(file, "%12.3f", complexes[num].coords[i].x); //31-38
		fprintf(file, "%8.3f", complexes[num].coords[i].y); //39-46
		fprintf(file, "%8.3f", complexes[num].coords[i].z); //47-54
		
		fprintf(file, "%6.2f", atom_occupancies[i]); //55-60
		fprintf(file, "%6.2f", atom_temp_factors[i]); //61-66
		fprintf(file, "%12s\n", at_symbols); //77-78
	}

	delete[] at_name;
	delete[] res_name;
	delete[] chain_id;
	delete[] at_symbols;
	delete[] pre_chain_cstring;
	
	fclose(file);
	return 0;
}
