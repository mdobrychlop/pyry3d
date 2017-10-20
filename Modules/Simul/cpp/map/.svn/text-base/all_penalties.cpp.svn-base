#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <iomanip>  
#include <iostream>
#include <cmath>
#include <cassert>
#include <tuple>
#include <sstream>
#include "all_penalties.hpp"
#include "clashes_map.hpp"
#include "technical.hpp"
#include "../pdb_parser/pdb_output_parser.hpp"

#define CA_CLASH_RADIUS 25.0
#define C4_CLASH_RADIUS 36.0

/* function are equivalent to python, only all_restraints are omitted*/

std::thread crysol_thread;

template <typename T> 
std::string to_str(const T& t) {
	/* helper function - translates anything for string (is translation is available) */
	std::ostringstream os; 
	os<<t; 
	return os.str(); 
}

void assign_new_taken_indexes(int complex_index) {

	if (needs_crysol())
		crysol_thread = std::thread(run_crysol, complex_index);
	
	//~ complexes[complex_index].taken_mapcells = 0;
	for (unsigned int i = 0; i < complexes[complex_index].components.size(); i++) {
		if (complexes[complex_index].components[i].changed) {
			complexes[complex_index].components[i].counted_taken_mapcells = 0;
			complexes[complex_index].components[i].taken_density = 0;
			std::vector<int> *indexes = &(complexes[complex_index].components[i].taken_mapcells);

			for (unsigned int j = 0; j < indexes->size(); j++) {
				complexes[complex_index].Grid_taken_mapcells[(*indexes)[j]]--;
				if (!complexes[complex_index].Grid_taken_mapcells[(*indexes)[j]]) {
					complexes[complex_index].components[i].taken_density +=
						Grid.grid_cells_density[(*indexes)[j]];
					complexes[complex_index].components[i].counted_taken_mapcells++;
				}
			}
		}	
	}
}

void assign_scores(int complex_index, bool flag) {

	if (!flag)
		complexes[complex_index].clashes =
			-(complexes[complex_index].all_clashes * 100.0 / complexes[complex_index].pairs_atoms); 
				 
	complexes[complex_index].outbox =
		-((complexes[complex_index].all_outbox_atoms * 100.0) / complexes[complex_index].alfa_atoms);

	int outatoms_cells = complexes[complex_index].mapcells_nr - complexes[complex_index].all_filled_cells;

	int empty_mapcells = Grid.mapcells - complexes[complex_index].all_filled_cells;
	
	int all_cells = complexes[complex_index].mapcells_nr + Grid.mapcells;
	
	if (freespace_penalty != 0){
		complexes[complex_index].freespace = -((((outatoms_cells + empty_mapcells) * 100.0) / all_cells));       
	} else {
		complexes[complex_index].freespace = 0.0f;
	}
	
	if (Grid.density_sum != 0)
		complexes[complex_index].density =
			-(((Grid.density_sum - complexes[complex_index].taken_densities) * 100.0)
			/ Grid.density_sum);
	else
		complexes[complex_index].density = 0.0;
	 
}

void __calculate_clashes_score(int complex_index, int component_index1, int component_index2) {
	
	bool need_to_compute = ( params.required_clashes_penalty ||
			params.required_clashes_penalty_allatoms );

	if (clashes_penalty && need_to_compute) {
		int clashes_num = 0;

		std::pair<int, int> conv1 = convert_chain_index_to_atom_start_end(component_index1);
		int atoms_start1 = conv1.first;
		int atoms_end1 = conv1.second;
		
		std::pair<int, int> conv2 = convert_chain_index_to_atom_start_end(component_index2);
		int atoms_start2 = conv2.first;
		int atoms_end2 = conv2.second;

		if(params.representation.compare("sphere") == 0 ){

			float rad1 = complexes[complex_index].coords[atoms_start2].atom_radii;
			float rad2 = complexes[complex_index].coords[atoms_start1].atom_radii;
			float radius = (1.01 * rad1) + rad2;
			//to się nie zgadza
			float sphere_distance = sqrt(dist_sqr(complexes[complex_index].coords[atoms_start2],
					complexes[complex_index].coords[atoms_start1]));

			float sphere_distance_penalty = fabs((rad1 + rad2) - sphere_distance);

			if(sphere_distance < radius){
				clashes_num ++;
			}

			complexes[complex_index].components[component_index1].sphere_dist_penalties[component_index2] =
				complexes[complex_index].components[component_index2].sphere_dist_penalties[component_index1] =
				sphere_distance_penalty;

		}else{

			for (int j = atoms_start1; j < atoms_end1; ++j){
				if(params.required_clashes_penalty &&  !is_alfa[j])
					continue;
				clashes_num += calculate_neighbour_clashes(j, 
						complexes[complex_index].coords[j], atoms_start2, atoms_end2,
						complex_index);
			}
		}

		complexes[complex_index].clashes_nr[component_index1][component_index2] = clashes_num;
		complexes[complex_index].clashes_nr[component_index2][component_index1] = clashes_num;

	}
}

void __calculate_empty_mapspaces_score(int complex_index, int component_index) {

	if (freespace_penalty || density_penalty || shapedesc) {
		
		complexes[complex_index].components[component_index].taken_density = 0;
		complexes[complex_index].components[component_index].taken_mapcells.clear();
		complexes[complex_index].components[component_index].counted_taken_mapcells = 0;
		std::pair<int, int> conv = convert_chain_index_to_atom_start_end(component_index);
		int atoms_start_index = conv.first;
		int atoms_end_index = conv.second;
		
		for (int i = atoms_start_index; i < atoms_end_index; i++) {
			
			float at_radii = atom_radii[atom_elem_symbols[i]];
			float radius = (Grid.radius / 2.0);
			
			//radius *= sqr3;
			if(params.representation.compare("sphere") == 0 ){
				radius = complexes[complex_index].coords[i].atom_radii + radius * sqr3;
			} else if (at_radii > radius) {
				radius *= sqr3;
				radius += at_radii;
			} else{
				radius *= sqr3;
			}

			__calculate_fitness_statistics(complex_index, complexes[complex_index].coords[i], radius, component_index);
		}
		
		complexes[complex_index].components[component_index].outmap_cells =
			complexes[complex_index].components[component_index].mapcells_nr
			- complexes[complex_index].components[component_index].taken_mapcells.size();

		if (complexes[complex_index].components[component_index].outmap_cells
			== complexes[complex_index].components[component_index].mapcells_nr) {
				
			complexes[complex_index].components[component_index].is_outmap = true;
			std::cerr << "COMPONENT OUTMAP!!!! True "<< complexes[complex_index].components[component_index].outmap_cells << std::endl;
		}
		else
			complexes[complex_index].components[component_index].is_outmap = false;
	}

}

void __calculate_fitness_statistics(int complex_index, point p, float radius,
		int component_index) {
	
	/* algorithm findins cells taken by atom with centre in p 
	 * and calculates statistics,
	 * Algorithm searchs grid_cell points in a square (2*radius) x (2*radius).
	 * For each grid_cell in this square checks
	 * if distance between grid_cell coords and p are <= radius
	 * */   
	
	float radius_sqr = radius * radius;
	
	int multi = ceil(radius / Grid.radius);
	point aligned = align_to_grid_cell_coord(p);
	float moved = multi * Grid.radius;
	
	float start_x = aligned.x - moved + Grid.radius;
	float start_y = aligned.y - moved + Grid.radius;
	float start_z = aligned.z - moved + Grid.radius;

	float end_x = aligned.x + moved;
	float end_y = aligned.y + moved;
	float end_z = aligned.z + moved;
	
	/* Sets boundary conditions */
	start_x = std::max(start_x, Grid.xmin);
	start_y = std::max(start_y, Grid.ymin); 
	start_z = std::max(start_z, Grid.zmin);
	
	/* Sets boundary conditions */
	end_x = std::min(end_x, Grid.xmax);
	end_y = std::min(end_y, Grid.ymax);
	end_z = std::min(end_z, Grid.zmax);

	float offset = Grid.radius / 2.0;
	/* end moved because of floating precision */
	end_x += offset;
	end_y += offset;
	end_z += offset;

	for (float i = start_x; i <= end_x; i+=Grid.radius) {
		float sqr_dist_x = pow(i - p.x, 2);
		for (float j = start_y; j <= end_y; j+=Grid.radius) {
			float sqr_dist_y = pow(j - p.y, 2);
			float sum = sqr_dist_x + sqr_dist_y;
			for (float k = start_z; k <= end_z; k+=Grid.radius) {
				float sqr_dist_z = pow(k - p.z, 2);

				if(sum + sqr_dist_z <= radius_sqr) {
					
					point curr(i, j, k);
					//~ if(dist_sqr(p, curr) <= radius_sqr) {
					int indx = get_gridcell_index_from_point(curr);

					if (Grid.grid_cells_ismap[indx]) {
						
						if (!complexes[complex_index].Grid_taken_mapcells[indx]) {

							complexes[complex_index].components[component_index].taken_density
								+= Grid.grid_cells_density[indx];
							complexes[complex_index].components[component_index].counted_taken_mapcells++;
						}
						
						complexes[complex_index].Grid_taken_mapcells[indx]++;
						complexes[complex_index].components[component_index].taken_mapcells.push_back(indx);	
						
					}
				}
			}
		}
	}

}

void calculate_outbox_mapfill(int complex_index) {
	
	for (unsigned int i = 0; i < complexes[complex_index].components.size(); i++) {
		
		if (complexes[complex_index].components[i].changed) {
		
			complexes[complex_index].all_outbox_atoms -= complexes[complex_index].components[i].outbox_atoms;
			complexes[complex_index].all_filled_cells -= complexes[complex_index].components[i].counted_taken_mapcells;
			complexes[complex_index].taken_densities -= complexes[complex_index].components[i].taken_density;
			///calculate new score elements
		
			__calculate_outbox_score(complex_index, i);
			__calculate_empty_mapspaces_score(complex_index, i);
			
			///add new scores to whole complex score
			
			complexes[complex_index].all_outbox_atoms += complexes[complex_index].components[i].outbox_atoms;
			complexes[complex_index].all_filled_cells += complexes[complex_index].components[i].counted_taken_mapcells;
			complexes[complex_index].taken_densities += complexes[complex_index].components[i].taken_density;
			
		}
	}
	
}

void __calculate_outbox_score(int complex_index, int component_index) {

	if (outbox_penalty != 0) { 
	
		std::pair<int, int> conv = convert_chain_index_to_atom_start_end(component_index);
		
		int atoms_start_index = conv.first;
		int atoms_end_index = conv.second;
		
		int outbox_atoms = 0;
		
		for (int i = atoms_start_index; i < atoms_end_index; i++) {
			if (is_alfa_atom(atom_names[i])) {
				if (!is_inside(complex_index, complexes[complex_index].coords[i])) {
					outbox_atoms++;
				}
			}
		}
		
		if (outbox_atoms >= 0.6 * complexes[complex_index].components[component_index].alfa_atoms) {
			
			std::pair<int, int> conv = convert_chain_index_to_atom_start_end(component_index);
			std::cerr << "COMPONENT OUTSIDE THE SIMULATION BOX "
			<< chain_identifiers[conv.first] << " " << outbox_atoms << " " 
			<< complexes[complex_index].components[component_index].alfa_atoms << std::endl;
			
			complexes[complex_index].components[component_index].is_outbox = true;
		} else
			complexes[complex_index].components[component_index].is_outbox = false;
		
		complexes[complex_index].components[component_index].outbox_atoms = outbox_atoms;
	}
}	

void calculate_restraints_clashes(int complex_index, bool all=false) {
	
	for (unsigned int i = 0; i < complexes[complex_index].components.size(); i++) {
		for (unsigned int j = i + 1; j < complexes[complex_index].components.size(); j++) {
			if(complexes[complex_index].components[i].changed || complexes[complex_index].components[j].changed || all) {
				complexes[complex_index].all_clashes -= complexes[complex_index].clashes_nr[i][j];
				__calculate_clashes_score(complex_index, i, j);
				complexes[complex_index].all_clashes += complexes[complex_index].clashes_nr[i][j];
				
			}
		}
	}
	
	float all_clashes_dist = get_all_clashes_dist(complex_index);

	complexes[complex_index].restraints -= all_clashes_dist;

}

void calculate_score_for_one_component_complex(int complex_index) {
	
	calculate_outbox_mapfill(complex_index);
	assign_scores(complex_index, true);

	if (needs_crysol())
		crysol_thread.join();

	complexes[complex_index].score = restraints_penalty * complexes[complex_index].restraints
		+ clashes_penalty * complexes[complex_index].clashes
		+ freespace_penalty * complexes[complex_index].freespace
		+ outbox_penalty * complexes[complex_index].outbox
		+ density_penalty * complexes[complex_index].density
		+ symmetry_penalty * complexes[complex_index].symmetry
		+ chi2_penalty * complexes[complex_index].chi2
		+ rg_penalty * complexes[complex_index].rg;
	
	std::cout << "ONE component score " << 
			complexes[complex_index].score << " components\t" <<
			( complexes[complex_index].restraints * restraints_penalty ) << '\t' <<
			( complexes[complex_index].clashes    * clashes_penalty ) << '\t' <<
			( complexes[complex_index].freespace  * freespace_penalty ) << '\t' <<
			( complexes[complex_index].outbox     * outbox_penalty ) << '\t' <<
			( complexes[complex_index].density    * density_penalty ) << '\t' <<
			( complexes[complex_index].symmetry   * symmetry_penalty ) << '\t' <<
			( complexes[complex_index].chi2       * chi2_penalty ) << '\t' <<
			( complexes[complex_index].rg  		  * rg_penalty ) << '\n';
}


float get_all_clashes_dist(int index) {
	float result = 0;
	int components_number = complexes[index].components.size();
	for(int i = 0; i < components_number; ++i ){
		for(int j = 0; j < components_number; ++j ){
			result += complexes[index].components[i].sphere_dist_penalties[j];
		}
	}
	return result;
}

void calculate_simulation_score_for_complex(int complex_index) {     

		int is_outbox = 0;
		int is_outmap = 1;
		
		for (unsigned int i = 0; i < complexes[complex_index].components.size(); i++) {
			if (complexes[complex_index].components[i].is_outbox) {
				is_outbox++;
				complexes[complex_index].is_complex_outbox = true;
			}
			
			if (complexes[complex_index].components[i].is_outmap) {
				is_outmap++;
				complexes[complex_index].is_complex_outmap = true;
			}
		}

		calculate_restraints_clashes(complex_index);
		assign_scores(complex_index, false);

		if(needs_crysol())
			crysol_thread.join();

		complexes[complex_index].score = 
								clashes_penalty * complexes[complex_index].clashes
								+ freespace_penalty * complexes[complex_index].freespace
								+ outbox_penalty * complexes[complex_index].outbox
								+ density_penalty * complexes[complex_index].density
								+ restraints_penalty * complexes[complex_index].restraints
								+ symmetry_penalty * complexes[complex_index].symmetry
								+ chi2_penalty * complexes[complex_index].chi2
								+ rg_penalty * complexes[complex_index].rg;
	   
		if (complexes[complex_index].is_complex_outbox) {
			std::cout << "OUTBOX COMPLEX!!!" << std::endl;
			float outbox_whole = std::max(100 * is_outbox * complexes[complex_index].outbox,
					-100000.0f);
			complexes[complex_index].score += outbox_whole;
		}
		
		if (complexes[complex_index].is_complex_outmap) {
			std::cout << "OUTMAP COMPLEX!!!" << std::endl;
			float outmap_whole = std::max(100 * is_outbox * (
					complexes[complex_index].density + complexes[complex_index].freespace),
					-100000.0f);
			complexes[complex_index].score += outmap_whole;
		}

		std::cout<< complexes[complex_index].score << "\t" <<
				( complexes[complex_index].restraints * restraints_penalty ) << '\t' <<
				( complexes[complex_index].clashes    * clashes_penalty ) << '\t' <<
				( complexes[complex_index].freespace  * freespace_penalty ) << '\t' <<
				( complexes[complex_index].outbox     * outbox_penalty ) << '\t' <<
				( complexes[complex_index].density    * density_penalty ) << '\t' <<
				( complexes[complex_index].symmetry   * symmetry_penalty ) << '\t' <<
				( complexes[complex_index].chi2       * chi2_penalty ) << '\t' <<
				( complexes[complex_index].rg  		  * rg_penalty ) << '\n';
}

bool is_inside(int complex_index, const point& p) {

	return (Grid.sb_xmin <= p.x 
			&& p.x <= Grid.sb_xmax
			&& Grid.sb_ymin <= p.y
			&& p.y <= Grid.sb_ymax
			&& Grid.sb_zmin <= p.z
			&& p.z <= Grid.sb_zmax);

}

bool needs_crysol(){
	return (rg_penalty > 0.0001 || chi2_penalty > 0.0001) && (params.curve.length() > 0);
}

void run_crysol(int complex_index) {

	//name of .pdb file
	std::string outname = params.name_prefix + to_str<float>(complexes[complex_index].score) + "_" + 
					to_str<int>(step_c) + "_" + to_str<float>(complexes[complex_index].temp) + ".pdb";

	write_to_file(complex_index, outname, atoms_number );

	std::string crysol_command = params.crysol_path + " " + outname + " " + params.curve;

	FILE *in;
	char buff[512];

	if(!(in = popen(crysol_command.c_str(), "r"))){
		std::cout<<"Wrong path for crysol given!";
		exit(EXIT_FAILURE);
	}

	std::string lines[6], line_with_result;
	int index_in_buffer = 0;
	while(fgets(buff, sizeof(buff), in)!=NULL){
		//fgets reads line by line
		lines[index_in_buffer] =std::string( buff );
		index_in_buffer = (index_in_buffer + 1) % 6;
	}

	pclose(in);

	line_with_result = lines[index_in_buffer];
	
	std::size_t chi_position = line_with_result.find("Chi:");
	std::size_t rg_position = line_with_result.find("RGT:");

	std::string chi_string = line_with_result.substr(chi_position + 4);
	std::string rg_string = line_with_result.substr(rg_position + 4, 6);

	float chi_notprocessed = std::stof(chi_string);
	float rg_notprocessed = std::stof(rg_string);

	if( chi_notprocessed <= 1.0 )
		chi_notprocessed = 0.0;
	else
		chi_notprocessed = - chi_notprocessed;

	float rge = - abs(params.rg_val - rg_notprocessed);

	int rm_success = system("rm *.log *.fit");

	if(chi2_penalty != 0)
		complexes[complex_index].chi2 = chi_notprocessed;
	if(rg_penalty != 0)
		complexes[complex_index].rg = rge;
}