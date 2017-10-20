#ifndef PYRY3D_CPP_HPP  
#define PYRY3D_CPP_HPP  

#include <cmath>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "complex.hpp"
#include "mutation/mutation.hpp"

#define CELL_BUFFER 14
/* number of cells left from each side of model for spare place  
TO DO: change python version and minimize this constant
*/

#define TRANSLATION_TYPES 9
/* rotation, translation, rotation of all, etc. */

#define CURR_COMPONENT(comp) (residue_indexes[ chain_indexes [ (comp) ] ])
/* a macro for finding index of first atom in given component */


#define NEXT_COMPONENT(comp) ( ((comp) >= chain_counter - 1) ? (atoms_number) : \
		(residue_indexes[ chain_indexes [ (comp) + 1 ] ]) );
/* a macro for finding index of first atom in next component */

struct covalent_bond{

	/* these indexes are relative to component, thay point at residues, not atoms!
	   at1_index and at2_index should be the same as in config file */
	int at1_index, at2_index;
	std::vector<int> comp_index;

};

struct disorder{
	std::pair < int, int > range;
	std::string type;
	float radius, max_sphere_radius;

	disorder(int begin, int end, std::string type, float radius,
			float max_sphere_radius): range(begin, end), type(type), 
			radius(radius), max_sphere_radius(max_sphere_radius) {}

	disorder() {}

};

//SINGLETON
struct components_common{

	std::vector < std::string > allowed_transformation; //potrzebne?

	std::vector < char > rot_axis; /* can contain X, Y, Z*/

	std::vector < covalent_bond > covalent_bonds; 

	std::vector < disorder > disorders;

	std::map < std::string, std::pair < float, float > > rot_ranges;
	std::map < std::string, std::pair < float, float > > rot_ranges_def;
	float trans_ranges[3][2];
	float trans_ranges_def[3][2];

	std::string moltype;
	float component_mass = 0;
};

#ifndef SWIG

//SINGLETON
struct transformation_freq {
	/* example config SETTINGS
	ROTATION_FREQ 0.5           
	ROTATION_COV_FREQ 0.2       
	TRANSLATION_FREQ 0.3
	the rest are 0s
	VECTORS:	frequencies [0,5; 0,8; 0,8; 0,8; 0,8; 0,8; 1; 1; 1]
				functions	[always point at the same, proper functions]
	*/
	std::vector<float> frequencies;

	void (* const functions[TRANSLATION_TYPES])(int, int, int, float&, char&,
			std::vector<float>&, bool, int&, int[]) = {&rotate, &translate, &exchange, &simul_dd,
			&translate, &rotate, &rotate, &rotate, &exchange};
	const char* names[TRANSLATION_TYPES] = {"rotation","translation","exchange","simul_dd",
			"translation_all","rotation_all","rotation_covalent", "rotation_whole", "exchange_and_sample"};
};

#endif

float dist(const point& a, const point& b);

void add_covalent_bond(std::vector<int> chain_indexes, int at1, int at2 );

void add_disorder(int start_pos, int stop_pos, std::string type, float radius,
			float max_sphere_radius);

void add_logical(int b, int f, bool and_);

void add_point_distance(const std::string& f, int a1, 
		int a2, float dist, const std::string& rel, float weig, 
		const std::vector<float>& point_v, const std::string& f_a_n);

void add_relation_interaction(const std::string& f, const std::string& s, 
		const std::string& t, const::std::string& fo, int a1, int a2, int b1, 
		int b2, int c1, int c2, int d1, int d2, const std::string& rel,
		const std::string& f_a_n, const std::string& s_a_n, const std::string& t_a_n,
		const std::string& fo_a_n);

void add_residue_distance(const std::string& f, const std::string& s, int a1, int a2, int b1, 
		int b2, float dist, const std::string& rel, float weig, const std::string& f_a_n,
		const std::string& s_a_n);

void add_surface_access(const std::string& f, int a1, int a2, 
		float dist, const std::string& rel, float weig, const std::string& f_a_n);

void add_symmetry(std::string f, std::string s, int a1, int a2, int b1, int b2,
		const std::string& f_a_n, const std::string& s_a_n);

void add_to_save(int index);

void append_score(int c);

void calculate_mass_of_complex_and_componets(int component_index);

void create_component(std::vector< std::string > transforms,
			std::vector< float > centre, std::vector < std::vector < float > >
			trans_ranges, std::vector < std::vector < float > > trans_ranges_sum,
			bool limited, std::vector< std::string > rot_axis,
			std::map < std::string, std::pair < float, float > > rot_ranges,
			std::map < std::string, std::pair < float, float > > rot_ranges_sum,
			std::string mol, bool dis, int all_components);

int deepcopy_complex(int index);
void deepcopy_complex_to_best_score_complex(int index);
void deepcopy_complex_to_best_score_niter_complex(int index);

void exchange(int i, int j);

void fill_can_be_exchanged();

void free_cpp();

int get_atoms_number();

float get_best_complex_score();

const std::vector<complex>& get_complexes();

std::vector<std::vector<float>> get_coords(int complex_index);

float get_niter_score();

float get_reptemp(int index);

int get_restraints_size();

void get_simul_params(float cl0, float cl1, float re0, float re1,
		float fr0, float fr1, float ou0, float ou1, float de0, float de1,
		float sy0, float sy1, float chi20, float chi21,
		float rg0, float rg1, std::string curve, std::string crysol_path,
		std::string name_prefix, float rg_val, std::string representation,
		bool required_clashes_penalty_allatoms,
		bool required_clashes_penalty);

float get_simulation_clashes(int index);
float get_simulation_density(int index);
float get_simulation_freespace(int index);
float get_simulation_outbox(int index);
float get_simulation_restraints(int index);
float get_simulation_score(int index);

int get_step();

float get_temp(int index);

void get_transformation_frequencies(std::vector < float > args);

float get_weigth(const std::string& symb);

void inc_step();

void init_cpp(int atoms, float clashes, float restraints, float freespace,
		float outbox, float density, float symmetry, float chi2, float rg, int iteration_count,
		bool scal, int chains, std::map<std::string, float> mols,
		bool shape, std::string logdir);

bool is_in_scores(int index);

int number_of_components();

int number_of_to_save();

int perform_mutation(int index, int iteration = -1);

int pool_size();

void pop();

	
void populate_pool(int max_pool_size);

void roulette(float max_temperature, int max_pool_size);

void save_best_niter(int iteration, std::string name);

std::vector<float> save_complexes_to_save(bool replica, const std::string& name);

void send_arrays(std::vector< float > range_bounds,
		std::vector< std::vector < float > > values,
		std::vector< int > mov, std::vector<int> fr,
		std::vector< std::vector < float > > surf);

void send_atom_radii(std::map<std::string, float>, std::vector <float> sphere_radii );

void send_grid(int mapcells, float density_sum, float sb_xmin, float sb_xmax, float sb_ymin, float sb_ymax,
	float sb_zmin, float sb_zmax, float xmin, float xmax, float ymin, float ymax,
	float zmin, float zmax, float radius, float overlap, std::vector<int> grid_cells_index,
	std::vector<float> grid_cells_density, std::vector<bool> grid_cells_ismap,
	std::vector<float> grid_cells_x, std::vector<float> grid_cells_y, std::vector<float> grid_cells_z);

void set_false_is_complex_outbox(int complex_index);

void set_false_is_complex_outmap(int complex_index);

void set_pool_size(int size);

void set_reptemp(int index, float value);

void set_temp (int index, float value);

int simul_size();			
	
void sort_by_simulation_score();

void tournament(float max_pool_size);

void try_to_push_in_to_save(int index);

#endif
