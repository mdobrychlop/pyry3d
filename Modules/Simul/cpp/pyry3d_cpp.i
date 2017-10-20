%module pyry3d_cpp
#include "complex.hpp"
%{
#define SWIG_FILE_WITH_INIT
#include "pdb_parser/pdb_input_parser.hpp"
#include "pdb_parser/pdb_output_parser.hpp"
#include "restraints/restraints.hpp"
#include "restraints/symmetry.hpp"
#include "variables.hpp"
#include "map/all_penalties.hpp"
#include "map/send_first_complex.hpp"
#include "map/technical.hpp"
%}

%include "std_map.i"
%include "std_pair.i"
%include "std_string.i"
%include "std_vector.i"

%include "carrays.i"
%array_functions(float, FloatArr)

namespace std {
	%template(StringFloatVectorMap)		map < string, pair < float, float > >;
	%template(FloatVector)				vector < float >;
    %template(FloatVectorVector)		vector < vector < float > >;
    %template(IntVector)				vector < int >;
    %template(StringVector)				vector < string >;
    %template(ComplexVector)			vector < complex >;
    %template(BoolVector)			vector < bool >;
    %template(StringFloatMap)			map < string, float >;
    %template(IntVectorVector)			vector < vector < int > >;
}


/* see interfejs.txt for explanation (in polish) */

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
		const std::string& f_a_n, const std::string& s_a_n, 
		const std::string& t_a_n, const std::string& fo_a_n);

void add_residue_distance(const std::string& f, const std::string& s, int a1, 
		int a2, int b1, int b2, float dist, const std::string& rel, float weig,
		const std::string& f_a_n, const std::string& s_a_n);

void add_surface_access(const std::string& f, int a1, int a2, 
		float dist, const std::string& rel, float weig, const std::string& f_a_n);

void add_symmetry(std::string f, std::string s, int a1, int a2, int b1, int b2,
		const std::string& f_a_n, const std::string& s_a_n);

void add_to_save(int index);


void append_score(int c);

void assign_new_taken_indexes(int complex_index);

void calculate_mass_of_complex_and_componets(int component_index);

void calculate_restraints_score(int index);

void calculate_symmetry_score(int index);

void create_component(std::vector< std::string > transforms,
			std::vector< float > centre, std::vector < std::vector < float > >
			trans_ranges, std::vector < std::vector < float > > trans_ranges_sum,
			bool limited, std::vector< std::string > rot_axis,
			std::map < std::string, std::pair < float, float > > rot_ranges,
			std::map < std::string, std::pair < float, float > > rot_ranges_sum,
			std::string mol, bool dis, int all_components);

void calculate_outbox_mapfill(int complex_index);

void calculate_simulation_score_for_complex(int complex_index);

void calculate_score_for_one_component_complex(int complex_index);

int deepcopy_complex(int index);
void deepcopy_complex_to_best_score_complex(int index);
void deepcopy_complex_to_best_score_niter_complex(int index);

void exchange(int i, int j);

void fill_can_be_exchanged();

void free_cpp();

int get_atoms_number();

float get_best_complex_score();

const std::vector<complex>& get_complexes();

std::vector<std::vector<float> > get_coords(int complex_index);

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

void inc_step();

void init_cpp(int atoms, float clashes, float restraints, float freespace,
		float outbox, float density, float symmetry, float chi2, float rg, int iteration_count,
		bool scal, int chains, std::map<std::string, float> mols,
		bool shape, std::string logdir);

bool is_in_scores(int index);

int number_of_components();

int number_of_to_save();

int parse_pdb(const std::string& file_name);

int perform_mutation(int index, int iteration = -1);

int pool_size();

void pop();

void set_false_is_complex_outbox(int complex_index);

void set_false_is_complex_outmap(int complex_index);

void set_pool_size(int size);
	
void populate_pool(int max_pool_size);

void roulette(float max_temperature, int max_pool_size);

std::vector<float> save_complexes_to_save(bool replica, const std::string& name);

void save_best_niter(int iteration, std::string name);

void send_arrays(std::vector<float> range_bounds,
		std::vector< std::vector < float > > values,
		std::vector<int> mov, std::vector<int> fr,
		std::vector< std::vector < float > > surf);

void send_atom_radii(std::map<std::string, float> ATOM_RADII, std::vector<float> sphere_radii);

void send_grid(int mapcells, float density_sum, float sb_xmin, float sb_xmax, float sb_ymin, float sb_ymax,
	float sb_zmin, float sb_zmax, float xmin, float xmax, float ymin, float ymax,
	float zmin, float zmax, float radius, float overlap, std::vector<int> grid_cells_index,
	std::vector<float> grid_cells_density, std::vector<bool> grid_cells_ismap,
	std::vector<float> grid_cells_x, std::vector<float> grid_cells_y, std::vector<float> grid_cells_z);

void set_changed_component(int complex_index, const std::vector<int>& changes);
void set_clashes_nr(int complex_index, std::vector<int> clashes);
void set_complex_data(int complex_index, int alfa_atoms,
	float all_clashes, int all_filled_cells, int all_outbox_atoms,
	float all_restraints, float clashes, float density,
	bool is_complex_outbox, bool is_complex_outmap, int mapcells_nr,
	int pairs_atoms, float restraints, float simulation_score,
	float symmetry, float taken_densities, int taken_mapcells);
	
void set_components_data(int complex_index, std::vector<bool> is_outbox,
	std::vector<bool> is_outmap, std::vector<int> mapcells_nr,
	std::vector<int> outbox_atoms);
	
void set_reptemp(int index, float value);

void set_temp (int index, float value);			

int simul_size();
	
void sort_by_simulation_score();

void tournament(float max_pool_size);

void try_to_push_in_to_save(int index);

void update_restraints_score(int index);

int write_to_file(int num, const std::string& file_name, int atoms_number);
