#include <string> 
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map> 
#include "technical.hpp"
#include "../pdb_parser/pdb_string_parser.hpp"

point align_to_grid_cell_coord(point p) {
	/* aligns atoms coords to coords of grid_cell */
	
	int multix = floor( (p.x - Grid.xmin) / Grid.radius );
	int multiy = floor( (p.y - Grid.ymin) / Grid.radius );
	int multiz = floor( (p.z - Grid.zmin) / Grid.radius );

	point res;
	res.x = Grid.xmin + (multix * Grid.radius); 
	res.y = Grid.ymin + (multiy * Grid.radius); 
	res.z = Grid.zmin + (multiz * Grid.radius);

	return res;
}

bool check_point(point p) {
	/* check if atom is inside Grid */
	if (p.x >= Grid.xmin && p.y >= Grid.ymin &&
		p.z >= Grid.zmin && p.x <= Grid.xmax &&
		p.y <= Grid.ymax && p.z <= Grid.zmax)
		return true;
	return false;
}

std::pair<int, int> convert_chain_index_to_atom_start_end(int chain_index) {
	/* return start and end position
	 * of chain_index in common complexes arrays */
	
	int residue_start_index = chain_indexes[chain_index];
	int residue_end_index;
	int atoms_start_index = residue_indexes[residue_start_index];
	int atoms_end_index;
	
	if (chain_index + 1 == chain_counter)
		atoms_end_index = atoms_number;
	else {
		residue_end_index = chain_indexes[chain_index + 1];
		atoms_end_index = residue_indexes[residue_end_index];
	}
	
	return std::make_pair(atoms_start_index, atoms_end_index);
}

float dist_sqr(const point& __restrict__ p1, const point& __restrict__ p2) {
	return pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2) + pow(p1.z - p2.z, 2);
}

int get_alfa_atoms(int component_index) {
	/* calculates alfta atoms for component 
	 * alfa: C4', C4*, CA
	 * */
	
	int res = 0;
	std::pair<int, int> conv = convert_chain_index_to_atom_start_end(component_index);
	
	int atoms_start = conv.first;
	int atoms_end = conv.second;

	for (int i = atoms_start; i < atoms_end; i++) {
		if (is_alfa_atom(atom_names[i]))
			res++;
	}
	
	return res;
}

int get_gridcell_index_from_point(point p) {
	
	/*takes grid_cell coords and return index of grid_cells*/

	int len_y = (round((Grid.ymax - Grid.ymin)
		/ Grid.radius)) + 1;
	
	int len_z = (round((Grid.zmax - Grid.zmin)
		/ Grid.radius)) + 1;
	
	int x = round((p.x - Grid.xmin) / Grid.radius);
	int y = round((p.y - Grid.ymin) / Grid.radius);
	int z = round((p.z - Grid.zmin) / Grid.radius);
	
	unsigned int res = (((x * (len_y)) + y) * (len_z)) + z;
	
	if (res >= complexes[0].Grid_taken_mapcells.size()) {
		std::cout << "res: [" << p.x << ", " << p.y << ", " << p.z << "] -> " << len_y << " -> " << len_z << std::endl;
	}
	return (((x * (len_y)) + y) * (len_z)) + z;
}

bool is_alfa_atom(const std::string& at_name) {
	return (at_name == "CA" or at_name == "C4'" or at_name == "C4*");
}
