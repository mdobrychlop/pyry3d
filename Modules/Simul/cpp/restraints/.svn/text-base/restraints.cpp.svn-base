#include "../extern.hpp"
#include "symmetry.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>

void calculate_restraints_score(int index){
	float score = 0.0f;

	int end = complexes[index].already_visited_by_logical.size();

	for(int k = 0; k < end; ++k){
		if(complexes[index].already_visited_by_logical[k])
			continue;
		int temp_i = k;
		for(int j = k + 1; j < end; ++j){
			// Visit only those, which were not visited earlier
			if(*(restraints[temp_i]) == *(restraints[j])){
				int res = calculate_diff(index, temp_i, j);
				if(res == 1)
					//we move to the correct score
					temp_i = j;
				complexes[index].already_visited_by_logical[j] = true;
			}
		}
		complexes[index].already_visited_by_logical[k] = true;
	}

    for(unsigned int i = 0; i < restraints.size(); ++i){

  		if(restraints[i]->type != 's')
	        if(complexes[index].res_scores[i] != -1){
	        	score += pow(complexes[index].res_scores[i], 2);
	        }
    }
    complexes[index].restraints = -1 * sqrt(score);


}

void update_restraints_score(int index){
	update_symmetry_score(index);
	struct complex *comp = &(complexes[index]);
	int i = 0;
	for(int k = 0; k < complexes[index].already_visited_by_logical.size(); ++k)
		complexes[index].already_visited_by_logical[k] = false;

	for(std::vector<Restraint*>::iterator it = restraints.begin();
			it != restraints.end(); ++it){
		if( (*it)->type!='s' ){
			if( (*it)->type=='q' ){
				RelationDistance* rel_dist = reinterpret_cast<RelationDistance*> (*it);
				rel_dist->score(index, i, rel_dist->first, rel_dist->second,
						rel_dist->third, rel_dist->fourth);	
			} else 
				(*it)->score(index, i, (*it)->first, (*it)->second, -1, -1);
				
		}
		++i;
	}
	for(std::vector<Logical>::iterator it = logical.begin(); it != logical.end();
			++it)
		it->score(index);
}


Logical::Logical(int b, int e, bool and_): id(logical.size()), begin(b), end(e), 
		max(and_){}

void Logical::score(int complex_index){
	/* it throws -1 everywhere where scores should not be inserted into result !!! 
	 Error prone*/

	float max_v = -1; 
	float min_v = std::numeric_limits<float>::max();
	int max_i = -1, min_i = -1;

	for(int i = begin; i < end; ++i){
		if(complexes[complex_index].already_visited_by_logical[i])
			continue;
		int temp_i = i;
		for(int j = i + 1; j < end; ++j){
			// Visit only those, which were not visited earlier
			if(*(restraints[temp_i]) == *(restraints[j])){
				int res = calculate_diff(complex_index, temp_i, j);
				if(res == 1)
					//we move to the correct score
					temp_i = j;
				complexes[complex_index].already_visited_by_logical[j] = true;
			}
		}
		complexes[complex_index].already_visited_by_logical[i] = true;
	}

	for(int i = begin; i < end; ++i){
		if(complexes[complex_index].res_scores[i] == -1)
			//artificial value
			continue;
		if(fabs(complexes[complex_index].res_scores[i]) > max_v){
			max_v = complexes[complex_index].res_scores[i];
			max_i = i;
		}
		if(fabs(complexes[complex_index].res_scores[i]) < min_v ){
			min_v = complexes[complex_index].res_scores[i];
			min_i = i;
		}
	}
	for(int i = begin; i < end; ++i){
		if(max){
			if(i != max_i)
				complexes[complex_index].res_scores[i] = -1;
		}else{
			if(i != min_i)
				complexes[complex_index].res_scores[i] = -1;
		}
	}

}

PointDistance::PointDistance(const std::string& s,  
		int a1, int a2,	float dist,	const std::string& rel, float weig, 
		const std::vector<float>& point_v, const std::string& f_a_n): 
		RelationRestraint(a1, a2, -2, -2, dist, rel, weig , 'r', f_a_n, f_a_n) {
	for(int i = 0; i < chain_counter; ++i){
		if(chain_identifiers[ residue_indexes [ chain_indexes[i] ] ] == s)
			first = i;
	}
	point new_point(point_v);
	this->p = new_point;
}
		

void PointDistance::get_score(int complex_index, int index, int first, int second){
	float min_dif = std::numeric_limits<float>::max();
	for(int i = this->first_range[0]; i < this->first_range[1]; ++i){
		int i_atom = residue_to_index[CURR_COMPONENT(first) + i - 1];
		while( atom_names[i_atom] != this->first_atom_name)
			/* we want to get atom which represents the whole residue */
			++i_atom;
		float act_dif = dist(
				complexes[complex_index].coords[ i_atom ],
				this->p) - this->distance;
		if(act_dif - std::numeric_limits<float>::epsilon() < 0){
			complexes[complex_index].res_scores[index] = 0;
			return;
		}
		if(act_dif < min_dif)
			min_dif = act_dif;
	}
	assert(min_dif < std::numeric_limits<float>::max());
	complexes[complex_index].res_scores[index] = this->weight * min_dif;
}

void RelationDistance::get_pair_score(int complex_index, int f, int s,
				int f_range[], int s_range[], float& min_diff, float& max_diff,
				const std::string& f_name, const std::string& s_name){
	min_diff = std::numeric_limits<float>::max();
	max_diff = -1;
	for(int i = f_range[0]; i < f_range[1]; ++i){
		int i_atom = residue_to_index[CURR_COMPONENT(f) + i - 1];
		while( atom_names[i_atom] != f_name )
			++i_atom;

		for(int j = s_range[0]; j < s_range[1]; ++j){
			int j_atom = residue_to_index[CURR_COMPONENT(s) + j - 1];
			while( atom_names[j_atom] != s_name )
				++j_atom;
			float act_dif = dist(
					complexes[complex_index].coords[ i_atom ],
					complexes[complex_index].coords[ j_atom ]);

			if(act_dif - std::numeric_limits<float>::epsilon() < 0)
				min_diff = 0;
			if(act_dif < min_diff)
				min_diff = act_dif;
			if(act_dif > max_diff)
				max_diff = act_dif;
		}
	}

	assert(min_diff < std::numeric_limits<float>::max());
	assert(max_diff >= 0);

}

RelationDistance::RelationDistance(const std::string& f, const std::string& s,
		const std::string& t, const std::string& fo,  
		int a1, int a2, int b1, int b2, int c1, int c2, int d1, int d2,
		const std::string& rel, const std::string& f_a_n, const std::string& s_a_n,
		const std::string& t_a_n, const std::string& fo_a_n): 
		Restraint(a1, a2, b1, b2, 'q', f_a_n, s_a_n), third_range{c1, c2+1}, 
		fourth_range{d1, d2+1}, third_atom_name(f_a_n), fourth_atom_name(fo_a_n){

	if(rel == "<")
		relation = 1;
	else if(rel == "<=")
		relation = 2;
	else if(rel == ">=")
		relation = 3;
	else if(rel == ">")
		relation = 4;
	else
		relation = 0;

	for(int i = 0; i < chain_counter; ++i){
		if(chain_identifiers[ residue_indexes [ chain_indexes[i] ] ] == f)
			first = i;
		if(chain_identifiers[ residue_indexes [ chain_indexes[i] ] ] == s)
			second = i;
		if(chain_identifiers[ residue_indexes [ chain_indexes[i] ] ] == t)
			third = i;
		if(chain_identifiers[ residue_indexes [ chain_indexes[i] ] ] == fo)
			fourth = i;
	}
}

void RelationDistance::get_score(int complex_index, int index, int first, int second,
		int third, int fourth){
	float res[2][2];
	get_pair_score(complex_index, first, second, first_range, second_range, 
		res[0][0], res[0][1], this->first_atom_name, this->second_atom_name);
	get_pair_score(complex_index, third, fourth, third_range, fourth_range, 
		res[1][0], res[1][1], this->third_atom_name, this->fourth_atom_name);
	if( ( this->relation == 1 ) || (this->relation == 2) )
		complexes[complex_index].res_scores[index] = ( res[0][0] <= res[1][1] ) ?
				0 : 100;
	else if( ( this->relation == 3 ) || (this->relation == 4) )
		complexes[complex_index].res_scores[index] = ( res[0][1] >= res[0][1] ) ?
				0 : 100;
}

void RelationRestraint::get_score(int complex_index, int index, int first, int second,
	int third, int fourth){
	this->get_score(complex_index, index, first, second);
}

int calculate_diff(int complex_index, int i,
		int j){

	RelationRestraint* res1 = (RelationRestraint*)(restraints[i]);
	RelationRestraint* res2 = (RelationRestraint*)(restraints[j]);

	int rel1 = res1->relation;
	int rel2 = res2->relation;

	rel1 = simplify_relation(rel1);
	rel2 = simplify_relation(rel2);


	if(rel1 == rel2){
		if (res1->weight * complexes[complex_index].res_scores[i] 
						< res2->weight * complexes[complex_index].res_scores[j]){
			complexes[complex_index].res_scores[j] = 0;
			return -1;
		}
		else {
			complexes[complex_index].res_scores[i] = 0;
			return 1;
		}

	} else if(rel1 == -rel2){
		if(rel1 == -1){
			if (res1->weight * complexes[complex_index].res_scores[i] 
						< res2->weight * complexes[complex_index].res_scores[j]){
				complexes[complex_index].res_scores[j] = 0;
				return -1;
			} else {
				complexes[complex_index].res_scores[i] = 0;
				return 1;
			}
		}
		else if(rel2 == -1){
			if (res2->weight * complexes[complex_index].res_scores[j] 
						< res1->weight * complexes[complex_index].res_scores[i]){
				complexes[complex_index].res_scores[i] = 0;
				return 1;
			} else {
				complexes[complex_index].res_scores[j] = 0;
				return -1;
			}
		}
	}

	//TO DO? example <, =
	return 0;
}

int simplify_relation(int rel){
	switch(rel){
		case 0:
			return 0;
			break;
		case 1:
			return -1;
			break;
		case 2:
			return -1;
			break;
		default:
			return 1;
	}
}

RelationRestraint::RelationRestraint(int a1, int a2, int b1, int b2, float dist, 
		const std::string& rel,	float weig, char t, const std::string& f_a_n,
		const std::string& s_a_n): Restraint(a1, a2, b1, b2, t, f_a_n, s_a_n),
		distance(dist), weight(weig) {
	if(rel == "<")
		relation = 1;
	else if(rel == "<=")
		relation = 2;
	else if(rel == ">=")
		relation = 3;
	else if(rel == ">")
		relation = 4;
	else
		relation = 0;
}

ResidueDistance::ResidueDistance(const std::string& s, const std::string& f, 
		int a1, int a2, int b1, int b2,
		float dist, const std::string& rel,	float weig, const std::string& f_a_n, 
		const std::string& s_a_n): RelationRestraint(a1, 
		a2, b1, b2, dist, rel, weig, 'd', f_a_n, s_a_n) {
	for(int i = 0; i < chain_counter; ++i){
		if(chain_identifiers[ residue_indexes [chain_indexes[i] ] ] == s)
			first = i;
		if(chain_identifiers[ residue_indexes [chain_indexes[i] ] ] == f)
			second = i;
	}
		
}


void ResidueDistance::get_score(int complex_index, int index, int first, int second){
	float min_dif = std::numeric_limits<float>::max();

	for(int i = this->first_range[0]; i < this->first_range[1]; ++i){
		int i_atom = residue_to_index[CURR_COMPONENT(first) + i - 1];

		while( atom_names[i_atom] != this->first_atom_name){
			/* we want to get atom which represents the whole residue */
			++i_atom;
		}
		for(int j = this->second_range[0]; j < this->second_range[1]; ++j){
			int j_atom = residue_to_index[CURR_COMPONENT(second) + j - 1];
			while( atom_names[j_atom] != this->second_atom_name)
				++j_atom;
			float act_dif = dist(
					complexes[complex_index].coords[ i_atom ],
					complexes[complex_index].coords[ j_atom ])
					- this->distance;
			if( ( ( (this->relation == 1) || (this->relation == 2) ) && (act_dif < 0) ) || 
					( ( (this->relation == 3) || (this->relation == 4) ) && (act_dif > 0) ) )
				min_dif = 0;
			/*if(act_dif - std::numeric_limits<float>::epsilon() < 0){
				min_dif = 0;
				break;
			}*/
			if(fabs(act_dif) < fabs(min_dif) )
				min_dif = act_dif;
		}
		if(min_dif - std::numeric_limits<float>::epsilon() < 0)
			break;
	}
	assert(min_dif < std::numeric_limits<float>::max());

	complexes[complex_index].res_scores[index] = this->weight * min_dif;
}

Restraint::Restraint(){}

Restraint::Restraint(int a1, int a2, int b1, int b2, char t, 
		const std::string& f_a_n, const std::string& s_a_n): type(t), 
		first_range{a1, a2+1}, second_range{b1, b2+1}, first_atom_name(f_a_n),
		second_atom_name(s_a_n) {}


void Restraint::score(int complex_index, int index, int first, int second,
			int third, int fourth){
	std::vector<int> args;
	args.push_back(first);
	if(second > -1)
		args.push_back(second);
	if(third > -1)
		args.push_back(third);
	if(fourth > -1)
		args.push_back(fourth);

	this->get_score(complex_index,  index,  first, second, third, fourth);
}

void RelationRestraint::print(){
	std::cout<<"RES1\n";
	std::cout<<"\tType : "<<(this->type)<<std::endl;
	std::cout<<"\tFirst : "<<(this->first)<<std::endl;
	std::cout<<"\tSecond : "<<(this->second)<<std::endl;

	std::cout<<"\tFirst atom name : "<<(this->first_atom_name)<<std::endl;
	std::cout<<"\tSecond atom name: "<<(this->second_atom_name)<<std::endl;
	std::cout<<"\tFirst range"<<this->first_range[0]<<" "<<this->first_range[1]<<std::endl;
	std::cout<<"\tSecond range"<<this->second_range[0]<<" "<<this->second_range[1]<<std::endl;
	std::cout<<"\tRelation"<<this->relation<<std::endl;
}

bool Restraint::operator==(const Restraint &q){	
	if (this->type == q.type){
		if((this->first_range[0] == q.first_range[0]) &&
				(this->first_range[1] == q.first_range[1]) &&
				(this->first == q.first) && 
				(this->first_atom_name == q.first_atom_name) &&
				(this->second_atom_name == q.second_atom_name)){
			if(this->type == 'r'){
				return true;
			}else if(this->type == 'd'){
				return ( (this->second == q.second) && 
					(this->second_range[0] == q.second_range[0]) &&
					(this->second_range[1] == q.second_range[1]) );
			}
		}
	}
	return false;
}


SurfaceAccess::SurfaceAccess(const std::string& s,
		int a1, int a2,	float dist, const std::string& rel,	float weig, const
		std::string& f_a_n): RelationRestraint(a1, a2, -2, -2, dist, rel, weig,
		'r', f_a_n, f_a_n) {
	for(int i = 0; i < chain_counter; ++i){
		if(chain_identifiers[ residue_indexes [ chain_indexes[i] ] ] == s)
			first = i;
	}
		
}


void SurfaceAccess::get_score(int complex_index, int index, int first, int second){
	float min_dif = std::numeric_limits<float>::max();
	for(int i = this->first_range[0]; i < this->first_range[1]; ++i){
		int i_atom = residue_to_index[CURR_COMPONENT(first) + i - 1];
		while( atom_names[i_atom] != this->first_atom_name)
			/* we want to get atom which represents the whole residue */
			++i_atom;

		for(unsigned int j = 0; j < surface.size(); ++j){
			float act_dif = dist(
					complexes[complex_index].coords[ i_atom ],
					surface[j])	- this->distance;
			if(fabs(act_dif) - std::numeric_limits<float>::epsilon() < 0){
				complexes[complex_index].res_scores[index] = 0;
				return;
			}
			if(fabs(act_dif) < min_dif)
				min_dif = fabs(act_dif);
		}
	}
	assert(min_dif < std::numeric_limits<float>::max());
	complexes[complex_index].res_scores[index] = this->weight * min_dif;
}