#ifndef CLASHES_MAP_PYRY3D_HPP 
#define CLASHES_MAP_PYRY3D_HPP 

#include "../pyry3d_cpp.hpp"

int calculate_neighbour_clashes(int atom_index, const point& __restrict__ p,
		int start2, int end2, int complex_index);

int get_index(const point& p);
	
void fill_containers(int complex);

bool outbox(const point& p);

void refresh_components(int complex);

void send_axes(int x_cells, int y_cells, int z_cells, float rad);

#endif
