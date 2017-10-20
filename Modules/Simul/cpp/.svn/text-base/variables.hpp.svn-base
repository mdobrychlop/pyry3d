

#ifndef VARIABLES_PYRY3D
#define VARIABLES_PYRY3D

#include <map>
#include <random>
#include <string>
#include <vector>
#include "pyry3d_cpp.hpp"
#include "restraints/restraints.hpp"
#include "restraints/symmetry.hpp"
#include "map/clashes_container.hpp"
#include "logger/logger.hpp"

/* yes, everything is global. Remember that our priority is performance */
/* pointers are arrays */

 grid Grid;
 std::map<std::string, float> atom_radii;
 /* for translating atom name to radius */
 float complex_mass = 0;
 int residue_counter = 0;	// number of resiudes
 int chain_counter = 0;		// number of components (component == chain !!)
 int poolsize = 1;			// number of complexes inside simualation (without new ones)
 int step_c = 0;			// number of step
 int steps; 				// number of steps in whole simulation

 int atoms_number;			// number of atoms
 int* atom_ids;				
 float* atom_occupancies;	
 float* atom_temp_factors;
 std::string* atom_elem_symbols;
 int* residue_to_index;		// translates number of residue to index in our coords arrays

 std::string* residue_name;
 int* residue_numbers;
 int* residue_indexes; 		// keeps first indexes of atoms

 std::string* chain_identifiers;
 int* chain_indexes; 		// keeps first indexes of residues		

 std::vector<complex> complexes, complexes_to_save;

 std::vector< Cells > containers;

 struct complex best_score_complex, best_score_niter_complex;
 
 std::vector<int> movable, free_components;	// they don't differ much, but they do!
 std::vector<int> can_be_exchanged;			// additional array for exchange function - see mutation.cpp
 std::vector<int> with_disorders;
 
 std::ranlux24 m_twister;	
 	//actually, ranlux24 is not mersenne twister - we use subtract with carry generator

	//random devices which are used often
 	// they are global so we don't create them before every use
 std::uniform_int_distribution<int> uint_dist;
 std::uniform_int_distribution<int> bool_int;
 std::uniform_real_distribution<float> zero_one;
 std::uniform_real_distribution<float> minusone_one;

 std::vector<Restraint*> restraints;
 std::vector<Logical> logical;

 std::vector<point> surface;	//additional array for counting surface access penalty

 std::vector<float> scores;

//simulation input variables

 simul_params params;

 float clashes_penalty, restraints_penalty, freespace_penalty, outbox_penalty,
		density_penalty, symmetry_penalty, rg_penalty, chi2_penalty;

 bool shapedesc;

 std::vector < float > scale_range_bounds;
 std::vector < std::vector < float > > scale_values;

 std::vector < components_common > comp_comm;

 transformation_freq transform_functions; 

 bool scaling;

 std::map<std::string, float> molweights;
 
 std::vector<bool> is_alfa;

 std::string* atom_names;


//temporary ! variables - don't use

 std::vector<point> atom_coords;

 Logger logger;

#endif
