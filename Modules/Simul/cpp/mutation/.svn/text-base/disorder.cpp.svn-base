#include "../pyry3d_cpp.hpp"
#include "../extern.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <list>
#include <random>

#define DISORDER_TRIES 100
#define ONEMINUSEPS 0.95

void anchored_reposition(grapes& g, int i){
	if(i < 1)
		/* have you ever heard of assertions? */
		return;

	point *root_grape, *anchor_grape;
	root_grape = &(g.list[i-1]);
	anchor_grape = &(g.list[i]);


	g.translate(root_grape->reversed_point(), true);

	float xAngle = atan2(anchor_grape->y, anchor_grape->z);
	g.rotate(xAngle, 'X', true);
	float yAngle = -atan2(anchor_grape->x, anchor_grape->z);
	g.rotate(yAngle, 'Y', true);
	
}

bool check_chained_distance(grapes& g){
	assert(!isnan(g.anchor.x));	
    float distance = dist(g.master, g.anchor);
    return distance < ((g.grape_count-1) * 2.0 * g.radius);
}

bool check_if_collide(int complex_index, int component_index, const grapes& g,
			int j){
	struct complex *curr = &(complexes[complex_index]);
	int end = NEXT_COMPONENT(component_index);
	for(int i = CURR_COMPONENT(component_index);
			i < end; ++i){

		/* a for loop over the rest of component */
		if(i == g.global_index)
			i+=g.grape_count;
			
		if(dist(g.list[j], curr->coords[i]) < g.radius + 1.5f)
				return true;
	}
	return false;
}


std::string choose_modeling_method(std::string moltype, std::string type){

	if (type == "simulated_volume"){
        if( moltype == "dna")
            return "chain";
        else
            return "random";
    }else if (type == "internal")
    	return "anchored_chain";
    else if( (type == "cterm") || (type == "nterm") ){
        if (moltype == "protein")
        	return "random";
        else
        	return "chain";
    }
    return "";
}

void choose_random_anchor_position(int complex_index, int anchors[], 
		point anchor_points[], std::vector< disorder>::iterator it){

	for(int i = 0; i < 2; ++i){
		if(anchors[i] == -1)
			return;
		set_anchor_points(complex_index, anchors[i], anchor_points[i]);

		std::uniform_int_distribution<int> axis_dist(1,3);
		int axis = axis_dist(m_twister);

		std::uniform_int_distribution<int> z_o(0,1);
		float vec = (((float) z_o(m_twister)) - 0.5f ) * 4 * it->radius;

		switch(axis){
			case 1:
				anchor_points[i].x += vec;
				break;
			case 2:
				anchor_points[i].y += vec;
				break;
			case 3:
				anchor_points[i].z += vec;
				break;

		}

	}

}

bool collides(const grapes& g, int i){

	if (dist(g.master, g.list[i]) > (g.master_radius - g.radius))
		return true;
	for(int j = 0; j < i; ++j){
		if (dist(g.list[j], g.list[i]) < (2* g.radius)){
			return true;
		}
	}
	return false;
}

void copy_results(const grapes& g, int complex_index){

	for(int i = 0; i < g.list.size(); ++i ){
		complexes[complex_index].coords[ i + g.global_index ].x = g.list[i].x;
		complexes[complex_index].coords[ i + g.global_index ].y = g.list[i].y;
		complexes[complex_index].coords[ i + g.global_index ].z = g.list[i].z;
	}

}

void end_anchored_chain(grapes& g){
	while(!g.history.empty()){
		anchored_history operation = g.history.top();
		switch(operation.name){
			case 'T':
				g.translate(operation.p.reversed_point());
				break;
			case 'R':
				g.rotate(-operation.angle, operation.axis);
				break;
		};
		g.history.pop();
	}
}

point find_chain_point(const point &p, float r, float max_angle,
		float min_angle = 0.0f, bool force = false ){
	float a, teta;
	if (force or (max_angle == min_angle) ){
	    a = 0.0f;
	    teta = max_angle;
	}
	else{
		int cnt; //co to znaczy?
		a = M_PI * 0.2f * minusone_one(m_twister);
		cnt = 0;
		do{
	    	teta = max_angle * minusone_one(m_twister);
	    }while ( (teta < min_angle) and (cnt++ < 10) );
	}
	point w;
	w.x = p.x + r * sin(teta) * cos(a);
	w.y = p.y + r * sin(teta) * sin(a);
	w.z = p.z + r * cos(teta);	
	return w;
}

point find_sphere_point(const point& p, float r){
	float alfa = zero_one(m_twister) * 2.0f * M_PI;
    float u = minusone_one(m_twister);
    point w;
    w.x = p.x + r * ( cos(alfa) * sqrt(1 - pow(u,2)) );
    w.y = p.y + r * ( sin(alfa) * sqrt(1 - pow(u,2)) );
    w.z = p.z + r * u;
    return w;
}

void generate(grapes& g, int complex_index, int component_index){

	g.copy_point(g.master, 0);


	if(g.generation_method == "anchored_chain"){
		if(! check_chained_distance(g) ){
			g.generation_method = "chain";
            std::cout<<"Warning! Anchored chain distance too big, \
            		switching to regular chain\n";
        }else{
	        g.copy_point(g.anchor, 1);
	        anchored_reposition(g, 1);
			std::uniform_real_distribution<float> pi_dist(0, M_PI);
			float r = pi_dist(m_twister);
			g.rotate(r, 'Z', true);

        }
	}

	for(int i = 1; i < g.grape_count; ++i){
		int iter_count = 0;
		int n;
		float angle;
		if(g.generation_method=="anchored_chain")
			anchored_reposition(g, i-1);
		while (iter_count++ < DISORDER_TRIES){
			if(g.generation_method == "random"){
				std::uniform_int_distribution<int> gen_dist(0, i-1);
				n = gen_dist(m_twister);
			}else if(g.generation_method == "chain")
				n = i - 1;
			else if(g.generation_method == "anchored_chain"){
				if(i == 1)
					continue;
				n = i - 2;
			}else
				n = -1;
			assert(n > -1);
			assert((i != 1) || (g.generation_method != "anchored_chain"));

			if (g.generation_method == "anchored_chain"){
				angle = get_max_angle(g.list[n], g.list[n+1], i, g);
				g.copy_point( find_chain_point(g.list[n], 2.0f * g.radius, 0.75f*angle,
						0.75f * fmin(angle, M_PI/4.0f), i==(g.grape_count-1) ), i);
			}
			else 
				g.copy_point(find_sphere_point(g.list[n], 2.0f * g.radius), i);

			if (!collides(g, i)){
				if (g.detect_collisions){
					if (component_index > -1)
						if (!g.detect_collisions(complex_index, component_index,
								g, i))
							//tu cuś więcej w argumentach
							break;
				}else
					break;
			}
		}

		
	}



    if (g.generation_method == "anchored_chain")
        end_anchored_chain(g);

    copy_results(g, complex_index);

}

float get_max_angle(const point& root, const point& anchor, int i, const grapes& g){
	float distance = dist(root, anchor);
	float anchor_radius = g.radius * 2.0f * ONEMINUSEPS;
	float root_radius = (g.grape_count - i) * 2.0f * g.radius;

	if(distance > anchor_radius + root_radius)
		return 0.0f;
	if(distance < fabs( root_radius - anchor_radius))
		return M_PI;
	float a = (pow(root_radius,2) - pow(anchor_radius,2) + pow(distance,2))
			/ (2.0f * distance);
	float height = sqrt(pow(root_radius,2) + pow(a,2));
	float max_angle = atan2(height, a);
	if(distance < a)
		return M_PI - max_angle;
	return max_angle;
}

int random_number(int i){
	std::uniform_int_distribution<int> temp_dist(0,i);
	return temp_dist(m_twister);
}

void set_anchor_points(int complex_index, int res_index, point& p){

	p.x = complexes[complex_index].coords[res_index].x;
    p.y = complexes[complex_index].coords[res_index].y;
    p.z = complexes[complex_index].coords[res_index].z;
}

void set_center(int complex_index, int component_index){
	float sum_of_masses = 0.0f;
	point sum_of_coords(0, 0, 0 );
	point *index = &(complexes[complex_index].coords[0]);

	int end = NEXT_COMPONENT(component_index);
	for(int i = CURR_COMPONENT(component_index); i < end; ++i){
		float w = get_weigth(atom_elem_symbols[i]);
		sum_of_masses += w;
		sum_of_coords += (index[i] * w);
	}
	complexes[complex_index].components[component_index].centre_of_mass[0] =
			sum_of_coords.x / sum_of_masses;
	complexes[complex_index].components[component_index].centre_of_mass[1] =
			sum_of_coords.y / sum_of_masses;
	complexes[complex_index].components[component_index].centre_of_mass[2] =
			sum_of_coords.z / sum_of_masses;

}

void set_volume_simulation_parameters(int complex_index, int component_index,
		grapes& g, std::vector< disorder>::iterator it){

	int anchors[2] = {-1,-1};
	point anchor_points[2];

	g.grape_count = it->range.second - it->range.first + 1;
	g.master_radius = it->max_sphere_radius;
	g.radius = it->radius;

	g.generation_method = 
			choose_modeling_method(comp_comm[component_index].moltype, it->type);

	if( it->type != "simulated_volume"){
		g.detect_collisions = &check_if_collide;
		//sort ?
		if(it->type == "nterm")
			anchors[0] = residue_to_index[CURR_COMPONENT(component_index) +
					it->range.second + 1];
		else if(it->type == "cterm")
			anchors[0] = residue_to_index[CURR_COMPONENT(component_index) + 
					it->range.first - 1];
		else if(it->type == "internal"){
			anchors[0] = residue_to_index[CURR_COMPONENT(component_index) + 
					it->range.first - 1];
			anchors[1] = residue_to_index[CURR_COMPONENT(component_index) + 
					it->range.second + 1];
		}

		while( (atom_names[anchors[0]] != "CA") && (atom_names[anchors[0]] != "C4")
				&& (atom_names[anchors[0] ] != "C4'") && (atom_names[anchors[0]] != "C4*") )
			++anchors[0];

		if(anchors[1] != -1)	
			while( (atom_names[anchors[1]] != "CA") && (atom_names[anchors[1]] != "C4")
					&& (atom_names[anchors[1] ] != "C4'") && (atom_names[anchors[1]] != "C4*") )
				++anchors[1];

		choose_random_anchor_position(complex_index, anchors, anchor_points, it);
		g.anchor = anchor_points[1];
	}
	else{
		struct component *comp = &(complexes[complex_index].components[component_index]);
		anchor_points[0].x = comp->centre_of_mass[0];
		anchor_points[0].y = comp->centre_of_mass[1];
		anchor_points[0].z = comp->centre_of_mass[2];
	}
	g.master = anchor_points[0];
	std::cout<<"VOLUME SIMULATIONS IIIII\n";

}

bool simulate_disorder(int complex_index, int component_index){
      
    if ( comp_comm[component_index].disorders.empty() ) 
    	return false;
	std::vector< disorder > temp_disorders;

	std::uniform_int_distribution<long long int> set_dist(1, pow(2,
			comp_comm[component_index].disorders.size() ) - 1);

	long long int set_representation = set_dist(m_twister);

	std::list<int> to_include;

	int index = 0;
	while(set_representation){
		if(set_representation % 2)
			to_include.push_back(index);
		set_representation /= 2;
		++index;
	}


	for(std::list<int>::iterator it = to_include.begin(); it != to_include.end();
		++it)
		temp_disorders.push_back(comp_comm[component_index].disorders[*it]);
	

	for(std::vector< disorder>::iterator it = temp_disorders.begin();
			it != temp_disorders.end(); ++it){
		int tmp = residue_to_index[ CURR_COMPONENT(component_index) + it->range.first ];

		grapes g(complexes[complex_index].coords[ tmp ], tmp);


		set_volume_simulation_parameters(complex_index, component_index, g, it);

		generate(g, complex_index, component_index);


	}


	/* can be optimized, we don't have to count every atoms in, we can use 
	 * deltas */
	set_center(complex_index, component_index);

	return true;

}

void simul_dd(int complex_index, int component_index, int func_index,
		float& angle, char& axis, std::vector<float>& l_vector, bool p, 
		int& covbound_index, int atoms[]){
	/* simulate disorder */

	std::cout<<"@@Simulate DD - first choice "<<
			chain_identifiers[residue_indexes[chain_indexes[component_index]]]<<
			std::endl;
    
    if ( simulate_disorder(complex_index, component_index) )
		complexes[complex_index].components[component_index].changed = true;
	else{
		std::uniform_int_distribution<int> dis_dist(0, with_disorders.size() -1);
		int new_index = dis_dist(m_twister);

		std::cout<<"@@Simulate DD - second choice "<<
				chain_identifiers[residue_indexes[chain_indexes[new_index]]]<<
				std::endl;

		simulate_disorder(complex_index, new_index); 
		complexes[complex_index].components[new_index].changed = true;
	}

}
