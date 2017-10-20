
#include <cassert>
#include <cstring>
#include <iostream>
#include <vector>


#include "complex.hpp"


point::point(const std::vector<float>& f): x(f[0]), y(f[1]), z(f[2]), atom_radii(0){
	//f is 3-element vector of coordinates
	assert(f.size() == 3);
}

complex::complex(): score(0.0f), temp(0.0f) {
	alfa_atoms = 0;
	all_clashes = 0;
	all_clashes_dist = 0;
	all_filled_cells = 0;
	all_outbox_atoms = 0;
	all_restraints = 0;
	clashes = 0.0;
	density = 0.0;
	freespace = 0.0;
	is_complex_outbox = false;
	is_complex_outmap = false;
	mapcells_nr = 0;
	outbox = 0;
	pairs_atoms = 0;
	taken_densities = 0.0;
	taken_mapcells = 0;
	rg = 0.0;
	chi2 = 0.0;
}


const bool by_simulation_score_descending::operator() (complex const &a, 
		complex const &b) const { 
	//method for comparing complexes
	return a.score > b.score;
}

std::ostream& operator<<(std::ostream& os, const point& dt)
{
    os << "" << dt.x << ' ' << dt.y << ' ' << dt.z;
    return os;
}