
#include <algorithm>
#include <cassert>
#include <cfloat>	
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>

#include <map>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "extern.hpp"
#include "map/technical.hpp"
#include "map/clashes_container.hpp"
#include "pyry3d_cpp.hpp"
#include "pdb_parser/pdb_string_parser.hpp"
#include "pdb_parser/pdb_output_parser.hpp"
//~ #include "mutation/disorder.hpp"


//private - not visible from python

float get_weigth(const std::string& symb){
	if( molweights.find(symb) == molweights.end())
		return 0.0f;
	else
		return molweights.at(symb);
}

void calculate_mass_of_complex_and_componets(int complex_index) {
	
	complex_mass = 0.0;
	for (int i = 0; i < complexes[complex_index].components.size(); i++) {
		std::pair<int, int> conv = convert_chain_index_to_atom_start_end(i);
		float component_mass = 0.0;
		for (int j = conv.first; j < conv.second; j++) {
			float w = get_weigth(atom_elem_symbols[j]);
			component_mass += w;
		}
		complex_mass += component_mass;
		comp_comm[i].component_mass = component_mass;
	}
}

template <typename T> 
std::string to_str(const T& t) {
	/* helper function - translates anything for string (is translation is available) */
	std::ostringstream os; 
	os<<t; 
	return os.str(); 
}

void actualize_reduction(const std::vector<bool>& is_chosen){
	/* clears pool, changes pool_size, should be run after reduction in
		genetic algorithm */
	std::vector<complex> new_pool;
	std::vector<Cells> new_containers;
	for(unsigned int i = 0; i < complexes.size(); ++i)
		if(is_chosen[i]){
			new_pool.push_back(std::move(complexes[i]));
			new_containers.push_back(std::move(containers[i]));
		}
	complexes.clear();
	containers.clear();
	for(std::vector<complex>::iterator it = new_pool.begin(); it != new_pool.end();
			++it)
		complexes.push_back(std::move(*it));

	for(auto it = new_containers.begin(); it != new_containers.end();
			++it)
		containers.push_back(std::move(*it));

	/* note std::move - it is a new thing added in c++11, it's very helpful
		pushing back a std::move of a structure, doesn't copy the structure !!!
		it only adds reference! */
	set_pool_size( complexes.size() );

	/* we also actualize size of cells vector - remember that every complex has
		its own vector */
	//cells.resize( complexes.size() );
}

void deepcopy_to_symbol(int index, complex &symbol){
	symbol = complexes[index];	
}

float dist(const point& a, const point& b){
	/* distance between two points */
	return sqrt(pow(a.x-b.x,2) + pow(a.y-b.y,2) + pow(a.z-b.z,2));
}

std::vector<std::vector<float>> get_coords(int complex_index) {
	/* return coords of all atoms from complex */
	std::vector<std::vector<float>> res;
	for (unsigned int i = 0; i < complexes[complex_index].coords.size(); i++) {
		float x = complexes[complex_index].coords[i].x;
		float y = complexes[complex_index].coords[i].y;
		float z = complexes[complex_index].coords[i].z;
		std::vector<float> coord;
		coord.push_back(x);
		coord.push_back(y);
		coord.push_back(z);
		res.push_back(coord);
	}
	return res;
}

float scale_parameter(float max, bool only_positive = false){
	//equivalent to scale_parameter
	float where_are_we = (float)step_c / (float)steps;
	unsigned int range = 0;
	for ( ; (range < scale_range_bounds.size()) &&
		(where_are_we >= scale_range_bounds[range]) ; ++range ){
	}
	if(range)
		--range;
	std::uniform_real_distribution<float> real_dist( 
			scale_values[range][0], scale_values[range][1]);
	float scaled_param = real_dist(m_twister) * max;
	if ( !only_positive){
		std::uniform_int_distribution<int> int_dist(0,1);
		if ( int_dist(m_twister) )
			scaled_param *= -1.0;
	} 
	return scaled_param;	
}

float calculate_scaled_weight(float up, float down){
	//equivalent to __calcuateScaledWeight
	float range = up - down;
	float new_penalty = range < 0 ? down - scale_parameter(fabs(range), true)
			: down + scale_parameter(fabs(range), true);
	return new_penalty;
	
}

void scale_weights(){
	clashes_penalty = calculate_scaled_weight(params.clashes[0], params.clashes[1]);
	restraints_penalty = calculate_scaled_weight(params.restraints[0], params.restraints[1]);
	freespace_penalty = calculate_scaled_weight(params.freespace[0], params.freespace[1]);
	outbox_penalty = calculate_scaled_weight(params.outbox[0], params.outbox[1]);
	density_penalty = calculate_scaled_weight(params.density[0], params.density[1]);
	symmetry_penalty = calculate_scaled_weight(params.symmetry[0], params.symmetry[1]);	
	chi2_penalty = calculate_scaled_weight(params.chi2[0], params.chi2[1]);
	rg_penalty = calculate_scaled_weight(params.rg[0], params.rg[1]);
}

//public

void add_covalent_bond(std::vector<int> chain_indexes, int at1, int at2 ){
	/* adds covalent bond to last component in 0 complex !!! */

	int index = complexes[0].components.size() - 1;
	comp_comm[index].covalent_bonds.emplace_back();
	std::vector< covalent_bond >::reverse_iterator it = 
			comp_comm[index].covalent_bonds.rbegin();
	it->comp_index = chain_indexes;
	it->at1_index = at1;
	it->at2_index = at2;
} 

void add_disorder(int start_pos, int stop_pos, std::string type, float radius,
		float max_sphere_radius){
	/* see add_covalent_bond */
	int index = complexes[0].components.size() - 1;

	comp_comm[index].disorders.emplace_back(start_pos -1, stop_pos -1, type, radius, max_sphere_radius);
}

void add_logical(int b, int f, bool and_){
	/* if !and_ its Or relation */
	Logical l(b, f, and_);
	logical.push_back(l);
	complexes[0].rearranged.push_back(true);
}

void add_point_distance(const std::string& f, int a1, 
		int a2, float dist, const std::string& rel, float weig, 
		const std::vector<float>& point_v, const std::string& f_a_n){
	PointDistance *n = 
			new PointDistance(f, a1, a2, dist, rel, weig, point_v, f_a_n);
	restraints.push_back(n);
	complexes[0].res_scores.emplace_back(0.0);
	complexes[0].already_visited_by_logical.push_back(false);

}

void add_relation_interaction(const std::string& f, const std::string& s, 
		const std::string& t, const::std::string& fo, int a1, int a2, int b1, 
		int b2, int c1, int c2, int d1, int d2, const std::string& rel,
		const std::string& f_a_n,  const std::string& s_a_n, 
		const std::string& t_a_n,  const std::string& fo_a_n){
	RelationDistance *n = 
			new RelationDistance(f, s, t, fo, a1, a2, b1, b2, c1, c2, d1, d2, rel,
					f_a_n, s_a_n, t_a_n, fo_a_n);
	restraints.push_back(n);
	complexes[0].res_scores.emplace_back(0.0);
	complexes[0].already_visited_by_logical.push_back(false);
}

void add_residue_distance(const std::string& f, const std::string& s, int a1, 
		int a2, int b1,	int b2, float dist, const std::string& rel,
		float weig, const std::string& f_a_n, const std::string& s_a_n){
	ResidueDistance *n = 
			new ResidueDistance(f, s, a1, a2, b1, b2, dist, rel, weig, f_a_n,
					s_a_n);
	restraints.push_back(n);
	complexes[0].res_scores.emplace_back(0.0);
	complexes[0].already_visited_by_logical.push_back(false);
}

void add_surface_access(const std::string& f, int a1, int a2, 
		float dist, const std::string& rel, float weig, const std::string& f_a_n){
	SurfaceAccess *n = new SurfaceAccess(f, a1, a2, dist, rel, weig, f_a_n);
	restraints.push_back(n);
	complexes[0].res_scores.emplace_back(0.0);
	complexes[0].already_visited_by_logical.push_back(false);
}

void add_symmetry(std::string f, std::string s, int a1, int a2, int b1, int b2,
		const std::string& f_a_n, const std::string& s_a_n){
	Symmetry *n = new Symmetry(f, s, a1, a2, b1, b2, f_a_n, s_a_n);
	restraints.push_back(n);
	complexes[0].res_scores.emplace_back(0.0);
	complexes[0].already_visited_by_logical.push_back(false);
}

void add_to_save(int index){
	// equivalent to self.complexes_to_save.append(deepcopy(self.mPool[index]))
	complexes_to_save.emplace_back();
	std::vector<complex>::reverse_iterator it = complexes_to_save.rbegin();
	deepcopy_to_symbol(index, (*it));
}

void append_score(int c){
	scores.push_back(complexes[c].score);
}

void create_component(std::vector< std::string > transforms,
			std::vector< float > centre, std::vector < std::vector < float > >
			trans_ranges, std::vector < std::vector < float > > trans_ranges_sum,
			bool limited, std::vector< std::string > rot_axis,
			std::map < std::string, std::pair < float, float > > rot_ranges,
			std::map < std::string, std::pair < float, float > > rot_ranges_sum,
			std::string mol, bool dis, int all_components){


	/* emplace back is c++11 novelty. It creates new element at the end of 
		container. It call constructor with arguments given to this function.
		In this example it is default constructor with no parameters */
	complexes[0].components.emplace_back();

	if(dis)
		/* component has disordered element */
		with_disorders.push_back(complexes[0].components.size() - 1);


	std::vector<component>::reverse_iterator it = complexes[0].components.rbegin(); 
	it->centre_of_mass[0] = centre[0];
	it->centre_of_mass[1] = centre[1];
	it->centre_of_mass[2] = centre[2];
	it->moves_limited = limited;
	for (int i = 0; i < 3; ++i){
		it->trans_history_sum[i] = 0.0;
		it->rot_history_sum[i] = 0.0;
	}

	it->rot_history_sum[3]=0.0;

	it->rot_ranges_sum[0][0] = rot_ranges_sum["X"].first;
	it->rot_ranges_sum[1][0] = rot_ranges_sum["Y"].first;
	it->rot_ranges_sum[2][0] = rot_ranges_sum["Z"].first;
	it->rot_ranges_sum[3][0] = rot_ranges_sum["L"].first;
	it->rot_ranges_sum[0][1] = rot_ranges_sum["X"].second;
	it->rot_ranges_sum[1][1] = rot_ranges_sum["Y"].second;
	it->rot_ranges_sum[2][1] = rot_ranges_sum["Z"].second;
	it->rot_ranges_sum[3][1] = rot_ranges_sum["L"].second;


	std::vector<components_common>::iterator it2 = comp_comm.begin() 
			+ complexes[0].components.size() - 1;


	if(trans_ranges_sum.size())
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 2; ++j){
				it2->trans_ranges[i][j] = trans_ranges[i][j];
				it2->trans_ranges_def[i][j] = it->trans_ranges_sum[i][j] = trans_ranges_sum[i][j];
			}


	it2->moltype = mol;
	it2->allowed_transformation = transforms;
	it2->rot_ranges = rot_ranges;
	it2->rot_ranges_def = rot_ranges_sum;

	for(std::vector< std::string >::iterator it_rot = rot_axis.begin(); 
			it_rot != rot_axis.end(); ++it_rot){
		assert(it_rot->size());
		it2->rot_axis.push_back(it_rot->at(0));
	}

	complexes[0].components.back().sphere_dist_penalties.resize(all_components);
}

int deepcopy_complex(int index){ 
	/* equivalent to deepcopy(self.mPool[index]) (somehow, because it also 
	 * prepares for adding to pool) 

		Deepcopy_complex creates new complex
		but created complex is not included in the pool !! It is similiar to the
		behaviour in python version
	 */ 

	complexes.emplace_back();
	containers.emplace_back(containers[index]);
	std::vector<complex>::reverse_iterator it = complexes.rbegin();

	assert(complexes.size() == containers.size());


	deepcopy_to_symbol(index, (*it));	

	/* returns index */

	return complexes.size() - 1;
}

void deepcopy_complex_to_best_score_complex(int index){ 
	/* equivalent to best_score_complex = deepcopy(self.mPool[index]) */
	deepcopy_to_symbol(index, best_score_complex);
}
void deepcopy_complex_to_best_score_niter_complex(int index){ 
	/* equivalent to best_score_complex_niter = deepcopy(self.mPool[index]) */
	deepcopy_to_symbol(index, best_score_niter_complex);
}

void exchange(int i, int j){ 
	/* equivalent to lTempComplex = self.mPool[i]
					 self.mPool[i] = self.mPool[j]
					 self.mPool[j] = lTempComplex
					 */
	std::swap(complexes[i], complexes[j]);
	std::swap(containers[i], containers[j]);

}

void fill_can_be_exchanged(){
	for(std::vector<int>::iterator it = free_components.begin();
			it != free_components.end(); ++it)
		if(comp_comm[*it].covalent_bonds.size() == 0)
			can_be_exchanged.push_back(*it);
}

void free_cpp(){
	/* we don't delete arrays inside complexes, they are freed by destructors */

	for(std::vector<Restraint*>::iterator it = restraints.begin(); it != restraints.end();
		++it){
		delete *it;
	} 

	delete [] atom_ids;
	delete [] atom_names;
	delete [] atom_occupancies;
	delete [] atom_temp_factors;
	delete [] atom_elem_symbols;
	delete [] residue_name;
	delete [] residue_numbers;
	delete [] residue_indexes;
	delete [] chain_identifiers;
	delete [] chain_indexes;
	delete [] residue_to_index;
	
}

int get_atoms_number(){
	return atoms_number;
}

float get_best_complex_score(){
	return best_score_complex.score;
}

const std::vector<complex>& get_complexes() {
	return complexes;    
}

float get_niter_score(){
	return best_score_niter_complex.score;
}

float get_reptemp(int index){
	return complexes[index].reptemp;
}

int get_restraints_size(){
	return restraints.size();
}

float get_simulation_clashes(int index) {
	return complexes[index].clashes;
}

float get_simulation_density(int index) {
	return complexes[index].density;
}

float get_simulation_freespace(int index) {
	return complexes[index].freespace;	
}
float get_simulation_outbox(int index) {
	return complexes[index].outbox;	
}

float get_simulation_restraints(int index) {
	return complexes[index].restraints;	
}

float get_simulation_score(int index){
  	return complexes[index].score;
}

void get_simul_params(float cl0, float cl1, float re0, float re1,
		float fr0, float fr1, float ou0, float ou1, float de0, float de1,
		float sy0, float sy1, float chi20, float chi21,
		float rg0, float rg1, std::string curve, std::string crysol_path,
		std::string name_prefix, float rg_val, std::string representation,
		bool req_clashes_penalty_all,
		bool req_clashes_penalty) {

	params.clashes[0] = cl0;
	params.clashes[1] = cl1;
	params.restraints[0] = re0; 
	params.restraints[1] = re1;
	params.freespace[0] = fr0; 
	params.freespace[1] = fr1;
	params.outbox[0] = ou0; 
	params.outbox[1] = ou1; 
	params.density[0] = de0; 
	params.density[1] = de1;
	params.symmetry[0] = sy0; 
	params.symmetry[1] = sy1;
	params.chi2[0] = chi20;
	params.chi2[1] = chi21;
	params.rg[0] = rg0;
	params.rg[1] = rg1;
	params.curve = curve;
	params.required_clashes_penalty_allatoms = req_clashes_penalty_all; 
	params.required_clashes_penalty = req_clashes_penalty;
	params.crysol_path = crysol_path;
	params.name_prefix = name_prefix;
	params.rg_val = rg_val;
	params.representation = representation;
}

int get_step(){
	return step_c;
}

float get_temp(int index){
	return complexes[index].temp;
}

void get_transformation_frequencies(std::vector < float > args){
	float fr = 0;
	//ASSUMPTION : sum of frequencies = 1
	//see: pyry3d_cpp.hpp structs
	for(unsigned int i = 0; i < args.size(); ++i)
		transform_functions.frequencies.push_back(fr+=args[i]);
	if((fr<0.99f) ||(fr>1.01f) ) {
		std::cout << "Frequencies should sum up to 1. Please change your config." << std::endl;
		throw "Frequencies should sum up to 1. Please change your config.";
	}
}

void inc_step(){
	++step_c;
}

void init_cpp(int atoms, float clashes, float restraints, float freespace,
		float outbox, float density, float symmetry, float chi2, float rg, int iteration_count, 
		bool scal, int chains, std::map<std::string, float> mols,
		bool shape, std::string logdir) { 

	atoms_number = atoms;

	atom_ids = new int[atoms_number];
	atom_names = new std::string[atoms_number];
	atom_coords.resize(atoms_number);
	atom_occupancies = new float[atoms_number];
	atom_temp_factors = new float[atoms_number];
	atom_elem_symbols = new std::string[atoms_number];
	residue_name = new std::string[atoms_number];
	residue_numbers = new int[atoms_number];
	residue_indexes = new int[atoms_number];
	chain_identifiers = new std::string[atoms_number];
	chain_indexes = new int[atoms_number];
	residue_to_index = new int[atoms_number];
	
	clashes_penalty = clashes;
	restraints_penalty = restraints;
	freespace_penalty = freespace;
	outbox_penalty = outbox;
	density_penalty = density;
	symmetry_penalty = symmetry;
	chi2_penalty = chi2;
	rg_penalty = rg;
	shapedesc = shape;
	steps = iteration_count;

	comp_comm.resize(chains);

	scaling = scal;
	std::random_device rd;
	m_twister.seed(rd());

	molweights = mols;

	std::cout.width(10);
	std::cout.precision(10);
	//std::ios_base::sync_with_stdio(false);

	logger.set_log(logdir);

}

bool is_in_scores(int index){
	return find(scores.begin(), scores.end(), complexes[index].score) != scores.end();
}

int number_of_components(){
	return chain_counter;
}

int number_of_to_save(){
	//equivalent to len(self.complexes_to_save)
	return complexes_to_save.size();
}

int perform_mutation(int index, int iteration){

	int new_index = deepcopy_complex(index);

	if(iteration > -1)
		complexes[new_index].iteration_index = iteration;

	if (index == 0)		
		/* if only first complex scales weigths, we are sure that scaling will
		be made one time and in every iteration */	
		scale_weights();	

	int moved_index = movable[uint_dist(m_twister)];
	float rand_value = zero_one(m_twister);
	int i = 0;

	for(; rand_value >= transform_functions.frequencies[i]; ++i);
		/* this for iterates over possible mutations and chooses correct one
			it is very c++ style of writing for loops, please don't mind it */

	for(int j = 0; j < chain_counter; ++j)
		/* assumption: new complex is created from a counted one or this is
		   zero iteration */
		complexes[new_index].components[j].changed = false;

	if(iteration == 0)
		for(int j = 0; j < chain_counter; ++j)
			complexes[new_index].components[j].changed = true;

	mutate(i, moved_index, new_index);

	return new_index;
}

int pool_size(){
	return poolsize;
}

void pop(){
	/* deletes wrong complex (which is at the end) */
	complexes.pop_back();
	containers.pop_back();
}
	
void populate_pool(int max_pool_size){ 
	//equivalent to populatePool

	int initial_pool_size = ceil((float)max_pool_size * 0.25);
	std::cout<< "Populating with initial " << initial_pool_size << " clones\n";
	for(int i = 1; i < initial_pool_size; ++i){
		deepcopy_complex(0);
		//cells.emplace_back();
		//cells.back().resize( cells.front().size() );
	}
	set_pool_size(initial_pool_size);
	
}	

void roulette(float max_temperature, int max_pool_size){
	float sum = 0;

	std::vector<float> scores;

	scores.reserve( complexes.size() );

	for(std::vector<complex>::iterator it = complexes.begin(); it !=
			complexes.end(); ++it)
		scores.push_back(it->score);

	std::sort(scores.begin(), scores.end());
	float quantile = scores[scores.size() / 4];



	for(std::vector<complex>::iterator it = complexes.begin(); it !=
			complexes.end(); ++it)
		if (it->score > quantile)
			sum += ( it->score - quantile );



	float current_prob = 0.0f;
	int index = 0;
	std::vector< float > wheel;
	std::vector< bool > is_chosen;
	wheel.resize( complexes.size() );
	is_chosen.resize( complexes.size() );
	for(std::vector<complex>::iterator it = complexes.begin(); it !=
			complexes.end(); ++it){
		current_prob += ( it->score > quantile ?
			(it->score - quantile) / sum : 0 );
		wheel[index] = current_prob;
		is_chosen[index] = false;
		++index;
	}

	/* probabilities ought to sum up to 1. */

	int chosen = 0;
	float draw_result;
	while (chosen < max_pool_size){
		index = -1;
		draw_result = zero_one(m_twister);
		for(unsigned int i = 0; i < wheel.size(); ++i)
			if(wheel[i] > draw_result){
				index = i;
				break;
			}
		if( (index > -1) && (!is_chosen[index]) ){
			is_chosen[index] = true;
			++chosen;
		}
	}

	actualize_reduction(is_chosen);

}	

void save_best_niter(int iteration, std::string name){
	if( find(scores.begin(), scores.end(), best_score_niter_complex.score) 
				== scores.end() ){
		scores.push_back(best_score_niter_complex.score);
		std::string outname = name + to_str<float>(best_score_niter_complex.score) + "_" + 
					to_str<int>(iteration) + "_" +
					to_str<float>(best_score_niter_complex.temp) + ".pdb";
		int s = complexes.size();
		complexes.push_back(best_score_niter_complex);
		write_to_file(s, outname, atoms_number );
		complexes.pop_back();
	}

}

std::vector<float> save_complexes_to_save(bool replica, const std::string& name){
	std::vector<float> res;
	for( std::vector<complex>::iterator it = complexes_to_save.begin();
			it != complexes_to_save.end(); ++it)
		if( find(scores.begin(), scores.end(), it->score) 
				== scores.end() ){
			scores.push_back(it->score);
			std::string outname = name + to_str<float>(it->score) + "_" + 
					to_str<int>(it->iteration_index)+"_";
			if ( replica )
				outname = outname + to_str<float>(it->reptemp) + ".pdb";
			else
				outname = outname + to_str<float>(it->temp) + ".pdb";
			complexes[0] = std::move(*it);

			write_to_file(0, outname, atoms_number);
			
			//send data for logfile.write_file
			res.push_back(complexes[0].iteration_index);
			res.push_back(complexes[0].score);
			res.push_back(complexes[0].restraints);
			res.push_back(complexes[0].clashes);
			res.push_back(complexes[0].freespace);
			res.push_back(complexes[0].outbox);
			res.push_back(complexes[0].density);
		} else {
			res.push_back(complexes[0].iteration_index);
			res.push_back(complexes[0].score);
			res.push_back(-1);
			res.push_back(-1);
			res.push_back(-1);
			res.push_back(-1);
			res.push_back(-1);
			
		}
	
	return res;
		
}

void send_arrays(std::vector< float > range_bounds,
		std::vector< std::vector < float > > values,
		std::vector<int> mov, std::vector<int> fr,
		std::vector< std::vector < float > > surf){
	scale_range_bounds = range_bounds;
	scale_values = values;
	movable = mov;
	free_components = fr;
	for(std::vector< std::vector < float > >::iterator it = surf.begin();
			it != surf.end(); ++it)
		surface.emplace_back(*it);

	std::uniform_int_distribution<int> dist(0, movable.size()-1);
	std::uniform_int_distribution<int> dist_bool(0, 1);
	std::uniform_real_distribution<float> dist_f(0.0, 1.0);
	std::uniform_real_distribution<float> dist_min(-1.0, 1.0);
	zero_one = dist_f;
	uint_dist = dist;
	bool_int = dist_bool;
	minusone_one = dist_min;

}

void send_atom_radii(std::map<std::string, float> ATOM_RADII, std::vector<float> radii) {
	/* send atoms radiuses from python */
	if(params.representation.compare("sphere") == 0){
		for(int i = 0; i < atoms_number; ++i){
			complexes[0].coords[i].atom_radii = radii[i];
		}
		return;
	}
	atom_radii = ATOM_RADII;	
}

void send_grid(int mapcells, float density_sum, float sb_xmin, float sb_xmax, float sb_ymin, float sb_ymax,
	float sb_zmin, float sb_zmax, float xmin, float xmax, float ymin, float ymax,
	float zmin, float zmax, float radius, float overlap, std::vector<int> grid_cells_index,
	std::vector<float> grid_cells_density, std::vector<bool> grid_cells_ismap,
	std::vector<float> grid_cells_x, std::vector<float> grid_cells_y,
	std::vector<float> grid_cells_z) {
	/* send Grid and grid_cells from python
	 * sets: ismap, density
	 * taken_mapcells for first complex
	 * sb_xmin - boundaries of SimulBox
	 * xmin - boundaries of Grid
	 *   */
	
	Grid.mapcells = mapcells;
	Grid.density_sum = density_sum;
	Grid.radius = 2 * radius;
	
	Grid.sb_xmin = sb_xmin;
	Grid.sb_xmax = sb_xmax;
	Grid.sb_ymin = sb_ymin;
	Grid.sb_ymax = sb_ymax;
	Grid.sb_zmin = sb_zmin;
	Grid.sb_zmax = sb_zmax;
	
	Grid.xmin = FLT_MAX;
	Grid.xmax = -FLT_MAX;
	Grid.ymin = FLT_MAX;
	Grid.ymax = -FLT_MAX;
	Grid.zmin = FLT_MAX;
	Grid.zmax = -FLT_MAX;

	for (unsigned int i = 0; i < grid_cells_index.size(); i++) {
		//complex[0] already exists!
		complexes[0].Grid_taken_mapcells.push_back(0);
		Grid.grid_cells_density.push_back(grid_cells_density[i]);
		Grid.grid_cells_ismap.push_back(grid_cells_ismap[i]);
		
		Grid.xmin = std::min(Grid.xmin, grid_cells_x[i]);
		Grid.ymin = std::min(Grid.ymin, grid_cells_y[i]);
		Grid.zmin = std::min(Grid.zmin, grid_cells_z[i]);
		
		Grid.xmax = std::max(Grid.xmax, grid_cells_x[i]);
		Grid.ymax = std::max(Grid.ymax, grid_cells_y[i]);
		Grid.zmax = std::max(Grid.zmax, grid_cells_z[i]);
		
	}
}

void set_false_is_complex_outbox(int complex_index) {
	complexes[complex_index].is_complex_outbox = false;
}

void set_false_is_complex_outmap(int complex_index) {
	complexes[complex_index].is_complex_outmap = false;
}

void set_pool_size(int size){ 
	/* as not every complex in vector is in pool , there can be a few 
	 * temporary complexes. Poolsize != complexes.size()*/
	poolsize = size;
}

void set_reptemp(int index, float value){ 
	complexes[index].reptemp = value;
}
	
void set_temp (int index, float value){
	complexes[index].temp = value;
}

int simul_size(){
	return complexes.size();
}

struct CmpIndices {
	bool operator() (int const &a, 
			int const &b) const { 
		//method for comparing complexes
		return complexes[a].score > complexes[b].score;
	}
};
	
void sort_by_simulation_score(){ 
	/* sort uses std::swap, so it works fast enough */
	std::vector<int> complexes_indices(complexes.size(), 0);
	for (int i = 0 ; i != complexes.size() ; i++)
		complexes_indices[i] = i;

	std::sort(complexes_indices.begin(), complexes_indices.end(), CmpIndices());
	std::sort(complexes.begin(), complexes.end(), by_simulation_score_descending());

	std::vector<Cells> temp_containers;

	for (int i = 0 ; i != complexes.size() ; i++) {
		/* build new cells table 
			No copying done, as I use move semantics */
		temp_containers.push_back(std::move(containers[complexes_indices[i]]));
	}

	containers.clear();
	for (int i = 0 ; i != complexes.size() ; i++) {
		/* build new cells table 
			No copying done, as I use move semantics */
		containers.push_back(std::move(temp_containers[i]));
	}

}


void tournament(float max_pool_size){

	std::vector<bool> is_chosen;
	is_chosen.resize( complexes.size() );
	for(int chosen = 0; chosen < max_pool_size; ++chosen){
		std::uniform_int_distribution<int> complex_dist(0, complexes.size() - 1);
		int contender_1 = -1, contender_2 = -1;
		do{
			contender_1 = complex_dist(m_twister);
			if(is_chosen[contender_1])
				contender_1 = -1;
		}while(contender_1 == -1 );
		do{
			contender_2 = complex_dist(m_twister);
			if(is_chosen[contender_2])
				contender_2 = -1;
		}while(contender_2 == -1 );

		if(complexes[contender_1].score > complexes[contender_2].score)
			is_chosen[contender_1] = true;
		else
			is_chosen[contender_2] = true;
	}

	actualize_reduction(is_chosen);
}

void try_to_push_in_to_save(int index){
	/* equivalent to 
	 * self.complexes_to_save.sort( key=lambda Complex: Complex.simulation_score)
		for el in self.complexes_to_save:
		    if aComplex.simulation_score > el.simulation_score:
			self.complexes_to_save.remove(el)
			self.complexes_to_save.append(deepcopy(aComplex))
			break
			*/ 
	std::sort(complexes_to_save.begin(), complexes_to_save.end(),
			by_simulation_score_descending());
	if(index > -1){
		for(std::vector<complex>::reverse_iterator it = complexes_to_save.rbegin();
				it != complexes_to_save.rend(); ++it){
			if(complexes[index].score > it->score){
				deepcopy_to_symbol(index, (*it));
				break;
			}		
		}
	}
}


