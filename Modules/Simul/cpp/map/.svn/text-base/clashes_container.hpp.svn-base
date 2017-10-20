#ifndef PYRY3D_CONTAINER
#define PYRY3D_CONTAINER


//for point definition
#include "../complex.hpp"

class Container {
	public:
		Container();
		Container(unsigned int atoms_number, unsigned int atoms_in_block);
		bool put(unsigned int atom_number, unsigned int position);
		bool remove(unsigned int atom_number, unsigned int position);
		int give_block();
		void fill_copy(std::vector<int>& copy, unsigned int position);
		void put_copy(std::vector<int>& copy, unsigned int position);
		void free_block(unsigned int position);
		bool collides(const std::vector<point>& coords, int component_start, 
				int component_end, int position, int atom_index) const;
	private:
		std::vector<int> atoms_indices;

		/* gives index of first atom in block nr (index of this vector) */
		std::vector<int> block_start;
		std::vector<int> indices_start;

		/* if true - the block is empty, and can be used */
		std::vector<bool> pool;

		/* we cycle through the pool to find an empty cell */
		int pool_index;

		int number_of_blocks;
		int block_size;

};

class Cells {
	public:
		Cells();
		Cells(unsigned int atoms_number, unsigned int cells_number);
		void put(unsigned int ind, unsigned int atom_number);
		void remove(unsigned int ind, unsigned int atom_number);
		bool collides(const std::vector<point>& coords, int component_start, 
				int component_end, int block_index, int atom_index) const;

	private:
		std::vector<Container> containers;

		/* indices of cells at containers vector */
		std::vector<int> used_container_index;

		/* position in given by used_container_index container - number of block*/
		std::vector<int> positions;

		/* variables used during copying atoms indices between containers 
			- we don't want to create small vectors during consecutive transfers,
			so instead we will use one, large enough to store indices during 
			the largest possible transfer - I believe it can be achieved easier by
			.reserve method, but the result will be the same*/
		std::vector<int> copy;
		int copy_size;
};

#endif