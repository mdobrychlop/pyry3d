
#include "../extern.hpp"
#include "../map/all_penalties.hpp"
#include "../map/clashes_map.hpp"
#include "../pyry3d_cpp.hpp"

#include "../complex.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
#include <random>
#include <vector>

std::map<char, int> axes = { {'X', 0 } ,{'Y', 1} ,{'Z', 2} ,{'L', 3} };

inline float radians(float degrees){
	//degrees to radians
	return M_PI*degrees/180.0;
}

inline float degrees(float radians){
	//radians to degrees
	return radians*180.0/M_PI;
}

//private 

void build_rotation_matrix(int points[], float& angle, float rotation_matrix[][4],
		int complex_index){

	struct complex *c = &(complexes[complex_index]);
	float uvw[3] = {c->coords[points[1]].x - c->coords[points[0]].x, 
					c->coords[points[1]].y - c->coords[points[0]].y, 
					c->coords[points[1]].z - c->coords[points[0]].z};
	float abc[3] = {c->coords[points[0]].x, c->coords[points[0]].y, c->coords[points[0]].z};
	float uvw2[3];

	if(points[0] == points[1])
		throw "To define line you need two different points in 3D space. Please change coordinated for \
				 COVALENT_BOND definition in configuration file";

	for(int i = 0; i < 3; ++i)
		uvw2[i] = uvw[i] * uvw[i];
	
	float cosT = cos(angle);
	float sinT = sin(angle);
	float oneMinusCosT = 1-cosT;
	
	float l2 = uvw2[0] + uvw2[1] + uvw2[2];
	float l = sqrt(l2);

	rotation_matrix[0][0] = (uvw2[0] + ((uvw2[1] + uvw2[2]) * cosT))/l2;
	rotation_matrix[0][1] = ((uvw[0]*uvw[1] * oneMinusCosT - uvw[2]*l*sinT))/l2;
	rotation_matrix[0][2] = ((uvw[0]*uvw[2] * oneMinusCosT + uvw[1]*l*sinT))/l2;
	rotation_matrix[0][3] = (((abc[0]*(uvw2[1] + uvw2[2]) - uvw[0]*(abc[1]*uvw[1] 
		+ abc[2]*uvw[2])) * oneMinusCosT + (abc[1]*uvw[2] - abc[2]*uvw[1])*l*sinT))/l2;
	
	rotation_matrix[1][0] = (uvw[0]*uvw[1] * oneMinusCosT + uvw[2]*l*sinT)/l2;
	rotation_matrix[1][1] = (uvw2[1] + ((uvw2[0] + uvw2[2]) * cosT))/l2;
	rotation_matrix[1][2] = ((uvw[1]*uvw[2] * oneMinusCosT - uvw[0]*l*sinT))/l2;
	rotation_matrix[1][3] = (((abc[1]*(uvw2[0] + uvw2[2]) - uvw[1]*(abc[0]*uvw[0] 
		+ abc[2]*uvw[2])) * oneMinusCosT + (abc[2]*uvw[0] - abc[0]*uvw[2])*l*sinT))/l2;
	
	rotation_matrix[2][0] = (uvw[0]*uvw[2] * oneMinusCosT - uvw[1]*l*sinT)/l2;
	rotation_matrix[2][1] = (uvw[1]*uvw[2] * oneMinusCosT + uvw[0]*l*sinT)/l2;
	rotation_matrix[2][2] = (uvw2[2] + ((uvw2[0] + uvw2[1]) * cosT))/l2;
	rotation_matrix[2][3] = (((abc[2]*(uvw2[0] + uvw2[1]) - uvw[2]*(abc[0]*uvw[0] 
		+ abc[1]*uvw[1])) * oneMinusCosT + (abc[0]*uvw[1] - abc[1]*uvw[0])*l*sinT))/l2;
		
}

void check_rotation_limits(int complex_index, int component_index, char& axis,
		float& angle){

	struct component *comp = &(complexes[complex_index].components[component_index]);
	struct components_common *comm = &(comp_comm[component_index]);

	if(comm->rot_ranges_def[std::string(axis,1)].first 
			>= comp->rot_history_sum[axes[axis]])
		comp->rot_ranges_sum[axes[axis]][0] = 0;
	else if(comm->rot_ranges_def[std::string(axis,1)].second 
			<= comp->rot_history_sum[axes[axis]])
		comp->rot_ranges_sum[axes[axis]][1] = 0;
	else{
		comp->rot_ranges_sum[axes[axis]][0] =
				-( angle + fabs(comp->rot_ranges_sum[axes[axis]][0]) );
		comp->rot_ranges_sum[axes[axis]][1] = 
				comp->rot_ranges_sum[axes[axis]][1] - angle;
		assert(fabs(comp->rot_ranges_sum[axes[axis]][1]) + fabs(comp->rot_ranges_sum[axes[axis]][0]) 
				<= fabs(2* comm->rot_ranges_def[std::string(axis,1)].first));
	}
}

void choose_rotation_points(int inpoints[2], int complex_index, int component_index,
			int& covbound_index){
	/* inpoints are indexes in main atoms array */
	struct components_common *comm = &(comp_comm[component_index]);

	if (comm->covalent_bonds.size() == 0){
		inpoints[0] = residue_indexes[chain_indexes[component_index]];
		int n = NEXT_COMPONENT(component_index);
		inpoints[1] = n - 1;
		std::cout<<"first/last atom "<<inpoints[0]<<" "<<inpoints[1]<<std::endl;
		covbound_index = -1;
	}else{
		std::uniform_int_distribution<int> covalent_dist(0,comm->covalent_bonds.size() - 1);
		covbound_index = covalent_dist(m_twister);
		struct covalent_bond *cov = &(comm->covalent_bonds[covbound_index]);
		inpoints[0] = residue_indexes[chain_indexes[component_index] + cov->at1_index];	
		inpoints[1] = residue_indexes[chain_indexes[component_index] + cov->at2_index];
		std::cout<<"selected bond! "<<covbound_index<<std::endl;
		std::cout<<"selected atom "<<inpoints[0]<<" "<<inpoints[1]<<" "
				<<cov->at1_index<<" "<<cov->at2_index<<std::endl;
	}

}

void exchange_components(int first_free_index, int second_free_index, 
		int complex_index, bool p){
	float angle = NAN;
	char axis = 0;
	int covbound_index = -1;
	int atoms[2] = {-1, -1};

	std::vector<float> l_vector1; 
	std::vector<float> l_vector2;

	l_vector1.resize(3);
	l_vector2.resize(3);

	for(int i = 0; i < 3; ++i){
		l_vector1[i] = complexes[complex_index].components[first_free_index].centre_of_mass[i] - 
			complexes[complex_index].components[second_free_index].centre_of_mass[i];
		l_vector2[i] = -l_vector1[i];
	}

	translate(complex_index, first_free_index, 1,
		angle, axis, l_vector2, p, covbound_index, atoms);

	translate(complex_index, second_free_index, 1,
		angle, axis, l_vector1, p, covbound_index, atoms);


}

void rotate_around_covalent_bond(float &angle, int inpoints[2], int &covbound_index,
				int complex_index, int component_index){
	float rangle = radians(angle);
	float rot[3][4];	/* rotation matrix */
	struct complex* comp = &(complexes[complex_index]);

	if(inpoints[0] == -1)
		choose_rotation_points(inpoints, complex_index, component_index,
				covbound_index);

	build_rotation_matrix(inpoints, rangle, rot, complex_index);

	int end = NEXT_COMPONENT(component_index);

	for (int i = residue_indexes[ chain_indexes [ component_index ] ]; 
			i < end ; ++i ){
		float xyz[3] = {comp->coords[i].x, comp->coords[i].y, comp->coords[i].z};
		comp->coords[i].x = rot[0][0] * xyz[0] + rot[0][1] * xyz[1] + rot[0][2] * xyz[2] + rot[0][3];
		comp->coords[i].y = rot[1][0] * xyz[0] + rot[1][1] * xyz[1] + rot[1][2] * xyz[2] + rot[1][3];
		comp->coords[i].z = rot[2][0] * xyz[0] + rot[2][1] * xyz[1] + rot[2][2] * xyz[2] + rot[2][3];
	}

	std::cout<<"Rotate "<<chain_identifiers[inpoints[0]]<<" "<<angle<<" "<<
		inpoints[0]<<" "<<inpoints[1]<<" "<<covbound_index<<std::endl;
}

void rotate_component(float &angle, char& axis, int complex_index, 
		int component_index, int func_index, int old_comp_index = 0, bool move = false ){
	float rangle = radians(angle);
	float rot[3][3][3] = {1, 0, 0, 0, float(cos(rangle)), float(sin(rangle)), 0, //if axis is X
			float(-sin(rangle)), float(cos(rangle)),
				float(cos(rangle)), 0, float(-sin(rangle)), 0, 1, 0, 			 //if axis is Y
			float(sin(rangle)), 0, float(cos(rangle)),
				float(cos(rangle)), float(sin(rangle)), 0, 						 //if axis is Z
			float(-sin(rangle)), float(cos(rangle)), 0, 0, 0, 1};
	int index;
	float centre_of_mass[3];

	std::cout<<"ROTATE: "<<chain_identifiers[residue_indexes[chain_indexes[component_index]]]
			<<" "<<angle<<" "<<axis<<std::endl;

	if (axis == 'X'){
		index = 0;
	}else if (axis == 'Y'){ 
		index = 1;  
	}else if (axis == 'Z'){  
		index = 2;
	}else 
		throw "No such axis !!";

	int end = NEXT_COMPONENT(component_index);

	float centre_of_complex_mass[3];

	if (move){

		for (int i = 0; i < 3; ++i)
			centre_of_mass[i] = complexes[complex_index].components[old_comp_index].centre_of_mass[i];

	} else {

		for (int i = 0; i < 3; ++i) {
			centre_of_complex_mass[i] = 0;
			for (int j = 0; j < complexes[complex_index].components.size(); ++j)
				centre_of_complex_mass[i] += comp_comm[j].component_mass *
						complexes[complex_index].components[j].centre_of_mass[i];
			centre_of_complex_mass[i] /= complex_mass;
		}

		for (int i = 0; i < 3; ++i)
			centre_of_mass[i] = func_index == 7 ? centre_of_complex_mass[i] :
					complexes[complex_index].components[component_index].centre_of_mass[i];
	}

	for (int i = residue_indexes[ chain_indexes [ component_index ] ]; i < end;
		++i){
		float temp_coords[3], temp_new_coords[3];
		temp_coords[0] = complexes[complex_index].coords[i].x - centre_of_mass[0];
		temp_coords[1] = complexes[complex_index].coords[i].y - centre_of_mass[1];
		temp_coords[2] = complexes[complex_index].coords[i].z - centre_of_mass[2];

		for (int j = 0; j < 3; ++j){
			temp_new_coords[j] = 0.0f;
			for (int k = 0; k < 3; ++k)
				temp_new_coords[j] += temp_coords[k] * rot[index][j][k];
		}
		complexes[complex_index].coords[i].x = temp_new_coords[0] + centre_of_mass[0];
		complexes[complex_index].coords[i].y = temp_new_coords[1] + centre_of_mass[1];
		complexes[complex_index].coords[i].z = temp_new_coords[2] + centre_of_mass[2];

	} 
}


void rotate_cov(int complex_index, int component_index, int func_index,
		float& angle, char& axis, std::vector<float>& l_vector, bool p,
		int& covbound_index, int atoms[], int old_comp_index){

	struct component *comp = &(complexes[complex_index].components[component_index]);

	comp->rot_history_sum[axes[axis]] += angle;

	if (func_index == 6){
		//covalent
		rotate_around_covalent_bond(angle, atoms, covbound_index, complex_index,
				component_index);
	}else{

		rotate_component(angle, axis, complex_index, component_index,
				func_index, old_comp_index, true);
		atoms = 0;
		covbound_index = -1;
	}
	comp->changed = true;

	if (comp->moves_limited)
		check_rotation_limits(complex_index, component_index, axis, angle);

}


void rotate_covalent_components(int component_index, int complex_index, 
	std::vector<float>& l_vector, bool p, float& angle, char& axis,
	int& covbound_index, int atoms[], int func_index){

	struct components_common *comm = &(comp_comm[component_index]);
	int i = -1;

	if(func_index==6){
		if (covbound_index==-1)
			return;
		for(std::vector<int>::iterator it = comm->covalent_bonds[covbound_index].comp_index.begin();
				it != comm->covalent_bonds[covbound_index].comp_index.end(); ++it){
			if(!complexes[complex_index].components[*it].changed)
				rotate_cov(complex_index, *it, func_index, angle, axis, l_vector, p,
						i, atoms, component_index);
		}
	}else{
		for(std::vector<covalent_bond>::iterator itb = comm->covalent_bonds.begin();
				itb != comm->covalent_bonds.end(); ++itb)
			for(std::vector<int>::iterator it = itb->comp_index.begin();
					it != itb->comp_index.end(); ++it)
				if(!complexes[complex_index].components[*it].changed){
					rotate_cov(complex_index, *it, func_index, angle, axis, l_vector, p,
							i, atoms, component_index);
				}
	}

}

float scale_parameter_range(float min, float max){
	//equivalent to scale_parameter_range
	float where_are_we = (float)step_c / (float)steps;
	unsigned int range = 0;
	for ( ; (range < scale_range_bounds.size()) &&
		(where_are_we >= scale_range_bounds[range]) ; ++range );
	
	if(range)
		--range;
	std::uniform_real_distribution<double> real_dist( 
			scale_values[range][0], scale_values[range][1]);
	//real_dist returns multiplier
	if (min > max){
			float temp = min;
			min = max;
			max = temp;
	}

	float middle = (max + min) / 2.0 ;
	float test = real_dist(m_twister);
	max -= middle;
	float scaled_param = test * max;
	std::uniform_int_distribution<int> int_dist(0,1);
	if ( int_dist(m_twister) )
		scaled_param *= -1.0; 
	return scaled_param + middle;
}

void scale_translation_vectors(int component_index, std::vector<float>& l_vector,
			int complex){
	
	float min_max[3][2];
	struct component *comp = &(complexes[complex].components[component_index]);
	struct components_common *comm = &(comp_comm[component_index]);

	min_max[0][0] = std::max(comp->trans_ranges_sum[0][0],
			comm->trans_ranges[0][0]);	
	min_max[0][1] = std::min(comp->trans_ranges_sum[0][1],
			comm->trans_ranges[0][1]);	
	
	min_max[1][0] = std::max(comp->trans_ranges_sum[1][0],
			comm->trans_ranges[1][0]);	
	min_max[1][1] = std::min(comp->trans_ranges_sum[1][1],
			comm->trans_ranges[1][1]);	
	
	min_max[2][0] = std::max(comp->trans_ranges_sum[2][0],
			comm->trans_ranges[2][0]);	
	min_max[2][1]= std::min(comp->trans_ranges_sum[2][1],
			comm->trans_ranges[2][1]);	

	l_vector.resize(3);

	if (scaling) {
		l_vector[0] = scale_parameter_range(min_max[0][0], min_max[0][1]);
		l_vector[1] = scale_parameter_range(min_max[1][0], min_max[1][1]);
		l_vector[2] = scale_parameter_range(min_max[2][0], min_max[2][1]);
	} else {
		std::uniform_real_distribution<float> dist_f1(min_max[0][0], min_max[0][1]),
				dist_f2(min_max[1][0], min_max[1][1]), 
				dist_f3(min_max[2][0], min_max[2][1]);
		l_vector[0] = dist_f1(m_twister);
		l_vector[1] = dist_f2(m_twister);
		l_vector[2] = dist_f3(m_twister);
	}
	
}

void translate_covalent_components(int component_index, int new_index, 
	std::vector<float>& l_vector, bool p, int& covbound_index){
	float f = NAN;
	char c = 0;
	int i = -1;
	struct components_common* comm = &(comp_comm[component_index]);

	for(std::vector<covalent_bond>::iterator itb = comm->covalent_bonds.begin();
			itb != comm->covalent_bonds.end(); ++itb)
		for(std::vector<int>::iterator it = itb->comp_index.begin();
				it != itb->comp_index.end(); ++it){
			if(!complexes[new_index].components[*it].changed)
				transform_functions.functions[1](new_index, *it, 1, f, c, l_vector,
						p, i, 0);
		}
}

//global

void exchange(int complex_index, int component_index, int func_index,
		float& angle, char& axis, std::vector<float>& l_vector, bool p, 
		int& covbound_index, int atoms[]){

	if(can_be_exchanged.size() < 2){
		std::cout<<"Wrong config given. You allowed exchanged mutation, while";
		std::cout<<"\nthere is no component pair which can be exchanged\n";
	}
	assert(can_be_exchanged.size()>1);

	struct component *comp = &(complexes[complex_index].components[component_index]),
		*comp2;
	int first_free_index;

	if((comp->moves_limited) || (comp_comm[component_index].covalent_bonds.size())){
		/* We have chosen a component which can't be moved */
		std::uniform_int_distribution<int> free_dist(0, can_be_exchanged.size()-1);
		first_free_index = can_be_exchanged[free_dist(m_twister)];
		comp = &(complexes[complex_index].components[first_free_index]);
	}else
		first_free_index = component_index;

	std::vector<int> movable_copy = can_be_exchanged;

	/* getting rid of first_free_index value, we don't want to sample the same
	 * component again*/
	movable_copy.erase(std::remove(movable_copy.begin(), movable_copy.end(), 
			first_free_index), movable_copy.end());

		
	std::uniform_int_distribution<int> free_dist(0, movable_copy.size()-1);
	int second_free_index = movable_copy[free_dist(m_twister)];
	comp2 = &(complexes[complex_index].components[second_free_index]);

	assert(first_free_index != second_free_index);

	exchange_components(first_free_index, second_free_index, complex_index,
			p);
		
	comp->changed = true;
	comp2->changed = true;

	int chosen_indices[2] = {first_free_index, second_free_index};

	if(func_index == 8)
		//exchange and sample
		sample(complex_index, chosen_indices);

}

void mutate(int func_index, int component_index, int new_index){
	float angle = NAN;
	char axis = 0;
	std::vector<float> l_vector;
	bool p = movable.size();

	int covbound_index = -1;
	int atoms[2] = {-1, -1};

	if ( (func_index == 4) || (func_index==5) || (func_index==7) ){
		//translate all or rotate all or rotate whole
		if (movable.size() != (unsigned int)chain_counter)
			return;

		for ( int comp = 0; comp < chain_counter; ++comp)
			transform_functions.functions[func_index](new_index, comp, 
					func_index, angle, axis, l_vector, p, covbound_index, atoms);

	} else if ( func_index == 1 ){
			//translate
		transform_functions.functions[func_index](new_index, component_index, 
					func_index, angle, axis, l_vector, p, covbound_index, atoms);
		translate_covalent_components(component_index, new_index, l_vector, p,
				covbound_index);

	} else if ( (func_index == 0) || (func_index==6) ){
			//rotate or covalent
		transform_functions.functions[func_index](new_index, component_index, 
				func_index, angle, axis, l_vector, p, covbound_index, atoms);
		rotate_covalent_components(component_index, new_index, l_vector, p,
				angle, axis, covbound_index, atoms, func_index );
	} else {
		transform_functions.functions[func_index](new_index, component_index, 
				func_index, angle, axis, l_vector, p, covbound_index, atoms);

	}

	/* refreshes clashes arrays*/
	refresh_components(new_index);

}

void rotate(int complex_index, int component_index, int func_index,
		float& angle, char& axis, std::vector<float>& l_vector, bool p, 
		int& covbound_index, int atoms[]){

	struct component *comp = &(complexes[complex_index].components[component_index]);
	struct components_common *comm = &(comp_comm[component_index]);

	if( (isnan(angle)) || (axis == 0) ){
		if (func_index == 6)
			axis = 'L';
		else{
			std::uniform_int_distribution<int> dist(0, comp_comm[component_index].rot_axis.size()-1);
			axis = comp_comm[component_index].rot_axis[dist(m_twister)];
		}
		float rotmin = std::max(comp->rot_ranges_sum[axes[axis]][0], comm->rot_ranges[std::string(1, axis)].first);
		float rotmax = std::min(comp->rot_ranges_sum[axes[axis]][1], comm->rot_ranges[std::string(1, axis)].second);
		if(scaling)
			angle = scale_parameter_range(rotmin, rotmax);
		else{
			std::uniform_real_distribution<float> dist_angle(rotmin, rotmax);
			angle = dist_angle(m_twister);
		}
	}

	comp->rot_history_sum[axes[axis]] += angle;

	if (func_index == 6){
		//covalent
		rotate_around_covalent_bond(angle, atoms, covbound_index, complex_index,
				component_index);
	}else{
		rotate_component(angle, axis, complex_index, component_index, func_index);
		atoms = 0;
		covbound_index = -1;
	}
	comp->changed = true;

	if (comp->moves_limited)
		check_rotation_limits(complex_index, component_index, axis, angle);

}

void sample(int complex_index, int chosen_indices[]){

	int current_iteration = step_c;

	float number_of_mutations = 0.01 * (steps - step_c);
	number_of_mutations = std::min(number_of_mutations, 100.0f);
	number_of_mutations = std::max(number_of_mutations, 10.0f);

	calculate_simulation_score_for_complex(complex_index);

	for(int i = 0; i < number_of_mutations; ++i){
		++current_iteration;
		std::cout<<"searching... "<<i<<std::endl;

		int index_in_array = bool_int(m_twister);

		//0 is rotate, 1 is translate - lookup the functions
		int translate = bool_int(m_twister);
		int new_index = deepcopy_complex(complex_index);

		mutate(translate, chosen_indices[index_in_array], new_index);

		complexes[new_index].components[chosen_indices[0]].changed = true;
		complexes[new_index].components[chosen_indices[1]].changed = true;

		//\\the function from Simul_cpp.py
		update_restraints_score(new_index);
		calculate_restraints_score(new_index);
		calculate_symmetry_score(new_index);
	
		if (number_of_components() == 1) {
		
			assign_new_taken_indexes(new_index);
			calculate_score_for_one_component_complex(new_index);
			
		} else {
			set_false_is_complex_outbox(new_index);
			set_false_is_complex_outmap(new_index);
			assign_new_taken_indexes(new_index);
			calculate_outbox_mapfill(new_index);
			calculate_simulation_score_for_complex(new_index);
		}

		if(complexes[complex_index].score < complexes[new_index].score)
			exchange(complex_index, new_index);

		//deletes complex with the worse cost
		pop();
	}
}

void translate(int complex_index, int component_index, int func_index,
		float& angle, char& axis, std::vector<float>& l_vector, bool p, 
		int& covbound_index, int atoms[]){
	if (l_vector.size() == 0)
		/* when size is 0, no translations of l_vector were made during this step,
		 * we want to assure every component uses the same vector */
		scale_translation_vectors(component_index, l_vector, complex_index);

	struct component *comp = &(complexes[complex_index].components[component_index]);
	struct components_common *comm = &(comp_comm[component_index]);

	std::cout<<"TRANSLATE: "<<chain_identifiers[residue_indexes[chain_indexes[component_index]]]<<" [ "
		<<l_vector[0]<<", "<<l_vector[1]<<", "<<l_vector[2]<<" ]"<<std::endl;


	for (int i = 0; i < 3; ++i) {
		comp->centre_of_mass[i] += l_vector[i];
	}


	/* here translation is made */
	int end = NEXT_COMPONENT(component_index);

	for (int i = residue_indexes[ chain_indexes [ component_index ] ]; i < end;
		++i){
		complexes[complex_index].coords[i].x += l_vector[0];		
		complexes[complex_index].coords[i].y += l_vector[1];
		complexes[complex_index].coords[i].z += l_vector[2];
	}
	
	comp->changed = true;

	if (comp->moves_limited){
		for(int i = 0; i < 3; ++i){
			comp->trans_history_sum[i] += l_vector[i];
			if(comm->trans_ranges_def[i][0] >= comp->trans_history_sum[i])
				comp->trans_ranges_sum[i][0] = 0;
			else if (comm->trans_ranges_def[i][1] <=  comp->trans_history_sum[i])
				comp->trans_ranges_sum[i][1] = 0;
			else{
				comp->trans_ranges_sum[i][1] -= l_vector[i];
				comp->trans_ranges_sum[i][0] = -(l_vector[i] +
						fabs(comp->trans_ranges_sum[i][0]));
				assert(fabs(comp->trans_ranges_sum[i][1]) +
						fabs(comp->trans_ranges_sum[i][0])
						<= fabs(2* comm->trans_ranges_def[i][0]) + 0.1);
				/* This should never happen while translating! */
			}
		}
	}

}
