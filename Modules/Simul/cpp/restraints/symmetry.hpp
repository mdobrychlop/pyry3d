#ifndef SYMMETRY_PYRY3D
#define SYMMETRY_PYRY3D

#include "restraints.hpp"

class Symmetry : public Restraint{
	public:

		void get_score(int complex_index, int index, int first, int second, int 
				third, int fourth);

		Symmetry();

		Symmetry(const std::string& s, const std::string& f, int a1, int a2, 
			int b1, int b2, const std::string& f_a_n, const std::string& s_a_n);
};

void calculate_symmetry_score(int index);

void update_symmetry_score(int index);

#endif