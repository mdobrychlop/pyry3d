#include "../extern.hpp"
#include "../pyry3d_cpp.hpp"

#include <cassert>
#include <limits>

void calculate_symmetry_score(int index){
	float product = 1.0f;
    float score = 0.0f;
    int sym = 0;
    for(unsigned int i = 0; i < restraints.size(); ++i)
        if(restraints[i]->type == 's'){
        	product *= complexes[index].res_scores[i];
        	++sym;
        }
    float average = pow(product, 1.0/(float)sym);
    for(unsigned int i = 0; i < restraints.size(); ++i)
    	if(restraints[i]->type == 's')
        	score += fabs(complexes[index].res_scores[i] - average);
    complexes[index].symmetry = -1.0f * score;
}

void update_symmetry_score(int index){
	struct complex *comp = &(complexes[index]);
	int i = 0;
	for(std::vector<Restraint*>::iterator it = restraints.begin();
			it != restraints.end(); ++it){
		if((*it)->type == 's')
			(*it)->get_score(index, i, (*it)->first, (*it)->second, -1, -1);
		++i;
	}
}

void Symmetry::get_score(int complex_index, int index, int first, int second,
		int third, int fourth){
	float min_dif = std::numeric_limits<float>::max();
	for(int i = this->first_range[0]; i < this->first_range[1]; ++i){
		int i_atom = residue_to_index[CURR_COMPONENT(first) + i - 1];
		while( atom_names[i_atom] != first_atom_name )
			++i_atom;

		for(int j = this->second_range[0]; j < this->second_range[1]; ++j){
			int j_atom = residue_to_index[CURR_COMPONENT(second) + j - 1];
			while( atom_names[j_atom] != second_atom_name )
				++j_atom;
			float act_dif = dist(
					complexes[complex_index].coords[ i_atom ],
					complexes[complex_index].coords[ j_atom ]);
			if(act_dif - std::numeric_limits<float>::epsilon() < 0){
				complexes[complex_index].sym_scores[index] = 0;
				return;
			}
			if(act_dif < min_dif)
				min_dif = act_dif;
		}
	}
	assert(min_dif < std::numeric_limits<float>::max());
	complexes[complex_index].res_scores[index] = min_dif;
}


Symmetry::Symmetry(){}

Symmetry::Symmetry(const std::string& s, const std::string& f, int a1, int a2, 
		int b1, int b2, const std::string& f_name , const std::string& s_name ):
		Restraint(a1, a2, b1, b2, 's', f_name, s_name){
	for(int i = 0; i < chain_counter; ++i){
		if(chain_identifiers[ residue_indexes [ chain_indexes[i] ] ] == s)
			first = i;
		if(chain_identifiers[ residue_indexes [ chain_indexes[i] ] ] == f)
			second = i;
	}
}