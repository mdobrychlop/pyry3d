#ifndef PYRY3D_ALL_PENALTIES_HPP  
#define PYRY3D_ALL_PENALTIES_HPP  

#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <thread>
#include <vector>
#include "../complex.hpp"

extern float clashes_penalty, restraints_penalty, freespace_penalty, outbox_penalty,
		density_penalty, symmetry_penalty;
extern std::vector<complex> complexes;
extern std::string* chain_identifiers; 
extern	std::string* atom_elem_symbols;
extern simul_params params;

void assign_new_taken_indexes(int complex_index);

void assign_scores(int complex_index, bool flag);

void __calculate_clashes_score(int complex_index);

void __calculate_empty_mapspaces_score(int complex_index, int component_index);

void __calculate_fitness_statistics(int complex_index, point p, float radius,
		int component_index);

void calculate_outbox_mapfill(int complex_index);

void __calculate_outbox_score(int complex_index, int component_index);

void calculate_restraints_clashes(int complex_index, bool all);

void calculate_score_for_one_component_complex(int complex_index);

//~ void calculate_simulation_score(int complex_index, bool shapedesc);

void calculate_simulation_score_for_complex(int complex_index);

bool is_inside(int complex_index, const point& p);

void run_crysol(int complex_index);

bool needs_crysol();

float get_all_clashes_dist(int complex);

#endif
