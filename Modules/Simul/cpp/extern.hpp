#include <random>
#include <string>
#include <vector>

#include "pyry3d_cpp.hpp"
#include "restraints/restraints.hpp"
#include "restraints/symmetry.hpp"
#include "map/clashes_container.hpp"
#include "logger/logger.hpp"

/* see variables.hpp for details */

 extern float complex_mass;

 extern std::map<std::string, float> atom_radii;
 extern grid Grid;
 
 extern int residue_counter ;
 extern int chain_counter ;
 extern int poolsize ;	
 extern int step_c;
 extern int steps; 

 extern int atoms_number;
 extern int* atom_ids;
 extern float* atom_occupancies;
 extern float* atom_temp_factors;
 extern std::string* atom_elem_symbols;
 extern int* residue_to_index;

 extern std::string* residue_name;
 extern int* residue_numbers;
 extern int* residue_indexes; // keep first indexes of atoms

 extern std::string* chain_identifiers;
 extern int* chain_indexes; // keep first indexes of residues		

 extern std::vector<complex> complexes, complexes_to_save;
 
 extern std::vector< Cells > containers;


 extern struct complex best_score_complex, best_score_niter_complex;
 
 extern std::vector<int> movable, free_components;
 extern std::vector<int> can_be_exchanged;
 extern std::vector<int> with_disorders;
 
 extern std::ranlux24 m_twister;
 extern std::uniform_int_distribution<int> uint_dist;
 extern std::uniform_int_distribution<int> bool_int;
 extern std::uniform_real_distribution<float> zero_one;
 extern std::uniform_real_distribution<float> minusone_one;

 extern std::vector<Restraint*> restraints;
 extern std::vector<Logical> logical;

 extern std::vector<point> surface;

 extern std::vector<float> scores;

 extern std::vector<bool> is_alfa;
 
//simulation input variables

 extern simul_params params;

 extern float clashes_penalty, restraints_penalty, freespace_penalty, outbox_penalty,
		density_penalty, symmetry_penalty, rg_penalty, chi2_penalty;
		
 extern bool 	shapedesc;
 
 extern std::vector<float> scale_range_bounds;
 extern std::vector< std::vector <float> > scale_values;

 extern std::vector < components_common > comp_comm;

 extern transformation_freq transform_functions; 

 extern bool scaling;

 extern std::map<std::string, float> molweights;

//temporary variables

 extern std::string* atom_names;

 extern std::vector<point> atom_coords;

 extern Logger logger;
