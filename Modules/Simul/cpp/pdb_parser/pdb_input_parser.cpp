#include <iostream>
#include <stdio.h>
#include <string>
#include <algorithm>
#include <cstring>
#include <sstream>
#include <vector>
#include "pdb_string_parser.hpp"

pdb_line::pdb_line(
	std::string record_name,
	std::string atom_ids,
	std::string atom_name,
	std::string residue_name,
	std::string chain_identifiers,
	std::string residue_numbers,
	std::string coordsx,
	std::string coordsy,
	std::string coordsz,
	std::string occupancy,
	std::string temperature_factor,
	std::string elem_symbol):

	record_name(record_name),
	atom_ids(std::stoi(atom_ids)),
	residue_name(residue_name),
	atom_name(atom_name),
	chain_identifiers(chain_identifiers),
	residue_numbers(std::stoi(residue_numbers)),
	coordsx(std::stof(coordsx)),
	coordsy(std::stof(coordsy)),
	coordsz(std::stof(coordsz)),
	occupancy(std::stof(occupancy)),
	temperature_factor(std::stof(temperature_factor)),
	elem_symbol(elem_symbol) {}

int parse_pdb(const std::string& pdb_string) {
	static std::vector<pdb_line> already_read;
	static size_t size = 0;
	static bool first_read = true;

	if (first_read) {
		already_read = parse_first_file(pdb_string);
		fill_all_array(already_read); // filling global arrs
		size = already_read.size();
	}

	std::vector<point> coords(size);
	
	if (first_read) {
		coords = atom_coords;
		first_read = false;
		fill_res_to_ind_array();
	}

	if (!size) {
		std::cerr<< "ERROR: already_read - vector with pdb_lines is empty\n";
		return -1;
	}

	std::vector< std::vector<int> > vec_clashes_nr;
	std::vector< std::vector<int> > vec_all_clashes_dist;
	vec_clashes_nr.resize(chain_counter);
	vec_all_clashes_dist.resize(chain_counter);
	for (int i = 0; i < chain_counter; i++) {
		vec_clashes_nr[i] = std::vector<int>(chain_counter, 0);
		vec_all_clashes_dist[i] = std::vector<int>(chain_counter, 0);
	}	
	
	complex created(coords, vec_clashes_nr, vec_all_clashes_dist);	
	complexes.push_back(created);
	
	return size;
}

void fill_all_array(const std::vector<pdb_line>& already_read) {

	std::string prev_chain = "AA"; // name should be shorter string
	int prev_residue = -1; // starting from 1
	size_t size = already_read.size();

	for (unsigned int i = 0; i < size; i++) {

		atom_ids[i] = already_read[i].atom_ids;

		atom_names[i] = already_read[i].atom_name;
		chain_identifiers[i] = already_read[i].chain_identifiers;
		atom_coords[i].x = already_read[i].coordsx;
		atom_coords[i].y = already_read[i].coordsy;
		atom_coords[i].z = already_read[i].coordsz;
		atom_occupancies[i] = already_read[i].occupancy;
		atom_temp_factors[i] = already_read[i].temperature_factor;
		atom_elem_symbols[i] = already_read[i].elem_symbol;
		residue_name[i] = already_read[i].residue_name;
		
		int curr_residue = already_read[i].residue_numbers;
		residue_numbers[i] = already_read[i].residue_numbers;

		std::string curr_chain = chain_identifiers[i];

		if ((prev_residue != curr_residue) ||
				(params.representation.compare("sphere") == 0) ||
				(curr_chain.compare(prev_chain))){

			residue_indexes[residue_counter] = i; // keep first indexes of atoms

			if ((curr_chain.compare(prev_chain)) ||
				(params.representation.compare("sphere") == 0)) {
				chain_indexes[chain_counter] = residue_counter; // keep first indexes of residues

				prev_chain = already_read[i].chain_identifiers;
				chain_counter++;
			}

			prev_residue = curr_residue;
			residue_counter++;
		}
	}
}

void fill_res_to_ind_array(){
	/* this array is bit complicated. If you want to know where is residue_number X
	 * from component Y, you just have to read residue_to_index[index_of_first_atom(Y) + 
	 * X -1 ]
	 * Please note that when I'm using it, I don't subtract one, because I prepared
	 * my residue_numbers earlier 
	 */
	int comp_i = 0, curr_comp = 0;
	int comp = NEXT_COMPONENT(comp_i);
	for(int i = 0; i < atoms_number; ++i)
		residue_to_index[i] = -1;

	for(int i = 0; i < atoms_number; ++i){

		if ( residue_to_index[residue_numbers[i] - 1 + curr_comp] == -1 )
			/* then it is first atom of this residue */
			residue_to_index[residue_numbers[i] - 1 + curr_comp] = i;
		if (i == comp - 1){
			++comp_i;
			curr_comp = comp;
			comp = NEXT_COMPONENT(comp_i);
		} 
	}
}
