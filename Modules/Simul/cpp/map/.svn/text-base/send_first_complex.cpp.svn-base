
//~ #include "../cuda/cuda.hpp"
#include "../extern.hpp"
#include "../pyry3d_cpp.hpp"
#include <algorithm>
#include <iostream>
#include <cmath>
#include <cstring>
#include <vector>
#include "clashes_map.hpp"
#include "technical.hpp"
#include "clashes_container.hpp"
#include "all_penalties.hpp"

#define CA_CLASH_RADIUS 25.0
#define C4_CLASH_RADIUS 36.0

void calculate_first_complex_taken_mapcells_stats(int complex_index, int i) {
	
	/* calculate statistics for first complex */
	
	complexes[complex_index].components[i].taken_density = 0;
	complexes[complex_index].components[i].counted_taken_mapcells = 0;
	complexes[complex_index].components[i].alfa_atoms = get_alfa_atoms(i);
	__calculate_empty_mapspaces_score(complex_index, i);
	complexes[complex_index].all_filled_cells += complexes[complex_index].components[i].counted_taken_mapcells;
	complexes[complex_index].taken_densities += complexes[complex_index].components[i].taken_density;
}

void create_cells(float max_rad){
	float x_cells, y_cells, z_cells;

	x_cells = ceil( (Grid.sb_xmax - Grid.sb_xmin) /
			max_rad ) + 2;
	y_cells = ceil( (Grid.sb_ymax - Grid.sb_ymin) /
			max_rad ) + 2;
	z_cells = ceil( (Grid.sb_zmax - Grid.sb_zmin) /
			max_rad ) + 2;
	//+2 as division is absolute

	long long int size = (x_cells + 2*CELL_BUFFER)*(y_cells + 2*CELL_BUFFER)*(z_cells + 2*CELL_BUFFER);
	//another +2 gives a spare cell from each side of model


	//assumption - signed number representation is two's complementary - -1 is equal to 255
	send_axes(x_cells, y_cells, z_cells, max_rad);


	containers.emplace_back(atoms_number, size);

	fill_containers(0);
}

void set_changed_component(int complex_index, const std::vector<int>& changes) {
	//by default components are not changed
	for(unsigned int i = 0; i < complexes[complex_index].components.size(); i++)
		complexes[complex_index].components[changes[i]].changed = false;
		
	for(unsigned int i = 0; i < changes.size(); i++) {
		//~ std::cout << "changes.size() " << changes.size() << " " << changes[i] << std::endl;
		complexes[complex_index].components[changes[i]].changed = true;
		//std::cout << "changed " << changes[i] << std::endl;	
	}
}

void set_complex_data(int complex_index, int alfa_atoms,
	float all_clashes, int all_filled_cells, int all_outbox_atoms,
	float all_restraints, float clashes, float density,
	bool is_complex_outbox, bool is_complex_outmap, int mapcells_nr,
	int pairs_atoms, float restraints, float simulation_score,
	float symmetry, float taken_densities, int taken_mapcells) {
	
	/* send data for first complex */
	
	complexes[complex_index].alfa_atoms = alfa_atoms;
	complexes[complex_index].all_clashes = all_clashes;
	complexes[complex_index].all_filled_cells = all_filled_cells;
	complexes[complex_index].all_outbox_atoms = all_outbox_atoms;
	complexes[complex_index].all_restraints = all_restraints;
	complexes[complex_index].clashes = clashes;
	complexes[complex_index].density = density;
	complexes[complex_index].is_complex_outbox = is_complex_outbox;
	complexes[complex_index].is_complex_outmap = is_complex_outmap;
	complexes[complex_index].mapcells_nr = mapcells_nr;
	complexes[complex_index].pairs_atoms = pairs_atoms;
	complexes[complex_index].restraints = restraints;
	complexes[complex_index].score = simulation_score;
	complexes[complex_index].symmetry = symmetry;
	complexes[complex_index].taken_densities = taken_densities;
}

void set_components_data(int complex_index, std::vector<bool> is_outbox, 
	std::vector<bool> is_outmap, std::vector<int> mapcells_nr,
	std::vector<int> outbox_atoms) {
	
	/* send data for first component */
	
	complexes[complex_index].all_filled_cells = 0;
	complexes[complex_index].taken_densities = 0;
		
	for (unsigned int i = 0; i < complexes[complex_index].components.size(); i++) {
		
		complexes[complex_index].components[i].is_outbox = is_outbox[i];
		complexes[complex_index].components[i].is_outmap = is_outmap[i];
			
		complexes[complex_index].components[i].mapcells_nr = mapcells_nr[i];
		complexes[complex_index].components[i].outbox_atoms = outbox_atoms[i];
		
		calculate_first_complex_taken_mapcells_stats(complex_index, i);
	}

	bool sphere = (params.representation.compare("sphere") == 0);

	for(int i = 0; i < atoms_number; ++i){
		is_alfa.push_back(is_alfa_atom(atom_names[i]));
		if(sphere){
			continue;
		}else if(params.required_clashes_penalty){
			if(is_alfa_atom(atom_names[i]))
				complexes[complex_index].coords[i].atom_radii = (atom_names[i] ==
						"CA")? CA_CLASH_RADIUS : C4_CLASH_RADIUS;
			else
				complexes[complex_index].coords[i].atom_radii = 0;
		}else if(params.required_clashes_penalty_allatoms)
			complexes[complex_index].coords[i].atom_radii = 
					pow(atom_radii[atom_elem_symbols[i]] * 2.01, 2);
	}

 	float max_rad;
 	if(params.required_clashes_penalty){
	 	for(int i = 0; i < atoms_number; ++i)
	 		if(atom_names[i] == "C4"){
	 			max_rad = 6.0f;
	 			break;
	 		}	
	 	max_rad = 5.0f;
	}
	else{
		max_rad = 0.0f;
		for(int i = 0; i < atoms_number; ++i)
			if(complexes[0].coords[i].atom_radii > max_rad)
				max_rad = complexes[0].coords[i].atom_radii;
		max_rad = sqrt(max_rad);
	}

	create_cells(max_rad);

	if(params.representation.compare("sphere") == 0 ){
		calculate_restraints_clashes(complex_index, true);
	}

}

void set_clashes_nr(int complex_index, std::vector<int> clashes) {

	for (unsigned int i = 0; i < clashes.size(); i+=3) {	
		int comp1_nr = clashes[i];
		int comp2_nr = clashes[i + 1];
		int clashes_nr = clashes[i + 2];
		complexes[complex_index].clashes_nr[comp1_nr][comp2_nr] = clashes_nr;
		complexes[complex_index].clashes_nr[comp2_nr][comp1_nr] = clashes_nr;
	}
	
}
