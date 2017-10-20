#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "clashes_container.hpp"

//for dist_sqrt declaration
#include "technical.hpp"



Container::Container() {}

Container::Container(unsigned int atoms_number, unsigned int atoms_in_block) {

	number_of_blocks = std::ceil(((float) atoms_number) / atoms_in_block);

	block_size = atoms_in_block;

	block_start.resize(number_of_blocks);
	indices_start.resize(number_of_blocks);
	pool.resize(number_of_blocks);

	/* on initialization, all cells are empty and available */
	std::fill(pool.begin(), pool.end(), true);

	for(int i = 0; i < block_start.size(); ++i){
		block_start[i] = i*block_size;
		indices_start[i] = i*block_size;
	}

	atoms_indices.resize(number_of_blocks * atoms_in_block);

	std::fill(atoms_indices.begin(), atoms_indices.end(), -1);

	pool_index = 0;
}

int Container::give_block() {
	/* should never run infinitely as we provided enough space */
	while(!pool[pool_index]){
		++pool_index;
		if(pool_index == number_of_blocks)
			pool_index = 0;
	}

	pool[pool_index] = false;
	return pool_index;
}

bool Container::put(unsigned int atom_number, unsigned int position) {
	int abs_position = block_start[position];

	assert (atoms_indices[indices_start[position]] == -1);

	atoms_indices[indices_start[position]] = atom_number;

	int a = indices_start[position];
	indices_start[position]++;

	assert(indices_start[position] == a+1);

	return indices_start[position] - abs_position + 1 == block_size; 
}

bool Container::remove(unsigned int atom_number, unsigned int position) {
	int abs_position = block_start[position];
	int after_atoms = indices_start[position];

	/* we don't check last as there can't be a situation when full block is 
		taken - we resize (move between conteners) in such situation */

	for(int i = abs_position; i < abs_position + block_size - 1; i++) {
		
		if(atoms_indices[i] == atom_number) {
			atoms_indices[i] = atoms_indices[after_atoms - 1];

			atoms_indices[after_atoms - 1] = -1;
			indices_start[position]--;

			if(block_size == 4)
				return indices_start[position] == block_start[position];

			return indices_start[position] - abs_position <= (block_size >> 2);
		} 

	}


	/* should never happen - it means there was no such atom in block! */
	assert(1==2);
	return true;


}

void Container::fill_copy(std::vector<int>& copy, unsigned int position) {
	unsigned int abs_position = block_start[position];
	auto it_from = atoms_indices.begin() + abs_position;
	std::copy(it_from, it_from + block_size, copy.begin()); 	
}


void Container::put_copy(std::vector<int>& copy, unsigned int position) {
	unsigned int abs_position = block_start[position];
	auto it_to = atoms_indices.begin() + abs_position;
	std::copy(copy.begin(), copy.begin() + block_size / 2, it_to);

	while(atoms_indices[abs_position] != -1)
		abs_position++;

	indices_start[position] = abs_position; 
}

void Container::free_block(unsigned int position) {
	unsigned int abs_position = block_start[position];
	unsigned int end = abs_position + block_size;

	for(; abs_position < end; ++abs_position) {
		if(atoms_indices[abs_position] == -1)
			break;
		atoms_indices[abs_position] = -1;		
	}
	pool[position] = true;
	indices_start[position] = block_start[position]; 
}

bool Container::collides(const std::vector<point>& coords, int component_start, 
				int component_end, int position, int atom_index) const {
	unsigned int abs_position = block_start[position];

	for(int i = abs_position; i < abs_position + block_size; ++i) {

		if (atoms_indices[i] == -1) {
			break;
		}
		if ( (atoms_indices[i] >= component_start) 
				&& (atoms_indices[i] < component_end) ) {

			if(dist_sqr(coords[atom_index], coords[atoms_indices[i]])
					< coords[atoms_indices[i]].atom_radii )
				return true;
		}
	}

	return false;
}

Cells::Cells() {}

Cells::Cells(unsigned int atoms_number, unsigned int cells_number ) {
	/* provides 4x needed space */
	/* we size down at 25% fill, so in the worst situation we should expect
		25%*block_size in every block. 
		Then atoms_number = 25%*block_size*number_of_blocks. The rest is 
		easy math */
	int spare_constant = 4;
	int block_size = 4;

	used_container_index.resize(cells_number);
	positions.resize(cells_number);

	for( int i = 0; i < std::ceil( std::log2(atoms_number) ) - 1; ++i ){
		containers.emplace_back(atoms_number * spare_constant, block_size);
		block_size *= 2;
	}

	/* -1 shows here incorrect value */
	std::fill(used_container_index.begin(), used_container_index.end(), -1);
	std::fill(positions.begin(), positions.end(), -1);

	/* block size is 2 times larger than atoms number, so the largest possible
		transfer is block_size / 4 */
	copy.resize(block_size / 4);

}

void Cells::put(unsigned int ind, unsigned int atom_number){
	int container_index = used_container_index[ind];

	if(container_index == -1){
		/* then we grab a block of size 4 */
		container_index = 0;
		used_container_index[ind] = 0;
		positions[ind] = containers[0].give_block();

	}

	Container* container = &containers[container_index];

	if(ind == 163582)
		std::cout<<ind<<" P "<<atom_number<<std::endl;

	if(container->put(atom_number, positions[ind])) {
		/* then we need to transfer */
		container->fill_copy(copy, positions[ind]);
		container->free_block(positions[ind]);
		++container_index;
		++used_container_index[ind];
		container = &containers[container_index];
		positions[ind] = container->give_block();
		container->put_copy(copy, positions[ind]);
	}

}

void Cells::remove(unsigned int ind, unsigned int atom_number) {
	int container_index = used_container_index[ind];
	Container* container = &containers[container_index];


	if(container->remove(atom_number, positions[ind])) {
		/* then we need to transfer */
		container->fill_copy(copy, positions[ind]);
		container->free_block(positions[ind]);
		--container_index;
		--used_container_index[ind];

		if(container_index == -1)
			return;

		container = &containers[container_index];
		positions[ind] = container->give_block();
		container->put_copy(copy, positions[ind]);
	}
}

bool Cells::collides(const std::vector<point>& coords, int component_start, 
				int component_end, int block_index, int atom_index) const {

	int used_container = used_container_index[block_index];
	int position = positions[block_index];

	if (used_container == -1)
		/* then this cell is empty */
		return false;

	return containers[used_container].collides(coords, component_start, 
			component_end, position, atom_index);

}
