
#include <cassert>
#include <vector>

#include "../pyry3d_cpp.hpp"
#include "technical.hpp"
#include "clashes_container.hpp"


#include <cstring>

int x,y,z;
float radius;


/*


Representation of cells in model

_________________________________________
|						/\
|						|							
|						|
|						|
|						| CELL_BUFFER of cells
|						|
|						|
|						|
|					   _\/________________
|<------------------->|
|					  |
|CELL_BUFFER of cells |		SIMULBOX
|					  |
|					  |
|
|



*/

bool outbox(const point& p);

int calculate_neighbour_clashes(int atom_index, const point& __restrict__ p,
		int start2, int end2, int complex_index){
	const std::vector<point>& coords = complexes[complex_index].coords;


	const Cells& cells_representation = containers[complex_index];
	
	/* this method has its own limitation - it doesn't count clashes outside 
		box */
	if(outbox(p))
		return 0;

	/* walks through 27 cells (26 neighbours and own) and searches for a clash */
	
	int start_cell = p.cell- (y * z + z + 1);
	for(int i = 0; i < 3; ++i){
		int curr_h = start_cell;
		for(int j = 0; j < 3; ++j){
			for(int k = curr_h; k < curr_h + 3; k += 1)
				if (cells_representation.collides(coords, start2, end2, k, atom_index) )
					return 2;

			curr_h += z;
		}
		start_cell += y*z;
	}
	return 0;
}

int get_index(const point& p){
	int x_i = floor( (p.x - Grid.sb_xmin) / radius);
	int y_i = floor( (p.y - Grid.sb_ymin) / radius);
	int z_i = floor( (p.z - Grid.sb_zmin) / radius);

	if(x_i < -CELL_BUFFER)
		x_i = -CELL_BUFFER;
	else if(x_i > x - 1 + CELL_BUFFER)
		x_i = x - 1;
	if(y_i < -CELL_BUFFER)
		y_i = -CELL_BUFFER;
	else if(y_i > y - 1 + CELL_BUFFER)
		y_i = y - 1;
	if(z_i < -CELL_BUFFER)
		z_i = -CELL_BUFFER;
	else if(z_i > z - 1 + CELL_BUFFER)
		z_i = z - 1;

	x_i += CELL_BUFFER;
	y_i += CELL_BUFFER;
	z_i += CELL_BUFFER;

	return  (x_i * y * z + y_i * z + z_i); 
}


void fill_containers(int complex){
	for(int i = 0; i < atoms_number; ++i){
		if( (params.required_clashes_penalty && (is_alfa[i])) ||
					(params.required_clashes_penalty_allatoms) ){
			int ind = get_index(complexes[complex].coords[i]);
			containers[complex].put(ind, i);
			complexes[complex].coords[i].cell = ind;
			if(ind <= 0){
				std::cout<<ind<<" OOPS\n";
				assert(2==3);

			}
			

		}
	}

}

bool outbox(const point& p){
	
	int x_i = floor( (p.x - Grid.sb_xmin) / radius);
	int y_i = floor( (p.y - Grid.sb_ymin) / radius);
	int z_i = floor( (p.z - Grid.sb_zmin) / radius);

	if(x_i < 0)
		return true;
	else if(x_i > x - 1)
		return true;
	if(y_i < 0)
		return true;
	else if(y_i > y - 1)
		return true;
	if(z_i < 0)
		return true;
	else if(z_i > z - 1)
		return true;

	return false;
}


void refresh_components(int complex){
	for(int i = 0; i < chain_counter; ++i)
		if(complexes[complex].components[i].changed){
			int end = NEXT_COMPONENT(i);
			int j = CURR_COMPONENT(i);
			for(; j < end; ++j)
				if( (params.required_clashes_penalty && (is_alfa[j])) ||
						(params.required_clashes_penalty_allatoms) ){
					int curr_cell = complexes[complex].coords[j].cell;
					int idx = get_index(complexes[complex].coords[j]);
					if(curr_cell != idx) {
						containers[complex].remove(curr_cell, j);
						containers[complex].put(idx, j);
						complexes[complex].coords[j].cell = idx;
					}
				}
		}
}



void send_axes(int x_cells, int y_cells, int z_cells, float rad){
	x = x_cells;
	y = y_cells;
	z = z_cells;
	radius = rad;
}
