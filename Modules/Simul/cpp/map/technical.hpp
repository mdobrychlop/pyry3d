#ifndef PYRY3D_TECHNICAL_HPP  
#define PYRY3D_TECHNICAL_HPP  

#include <fstream>
#include <utility>
#include <vector>
#include "../complex.hpp"
#include "../extern.hpp"

const float sqr3 = 1.7320;


point align_to_grid_cell_coord(point p);

void alloc_coords(float* coordx, float* coordy, float* coordz, int atoms_in_chain,  int start_atom_index, int complex_index);

void calculate_stats_for_corner(int complex_index, int component_index,
	float x, float y, float z, point p, float radius_sqr);

void calculate_stats_for_edge_x(int complex_index, int component_index,
	float x, float y, float z, float end_x, float end_y, float end_z,
	point p, float radius_sqr); // move x

void calculate_stats_for_edge_y(int complex_index, int component_index,
	float x, float y, float z, float end_x, float end_y, float end_z,
	point p, float radius_sqr); // move y

void calculate_stats_for_edge_z(int complex_index, int component_index,
	float x, float y, float z, float end_x, float end_y, float end_z,
	point p, float radius_sqr); // move z

void calculate_stats_for_index(int complex_index, int component_index, int indx);

void calculate_stats_for_outer_points(int complex_index, int component_index,
	point p, float radius_sqr);

void calculate_stats_for_surface_x(int complex_index, int component_index,
	float x, float y, float z, float end_x, float end_y,
	float end_z, point p, float radius_sqr); // move y, z

void calculate_stats_for_surface_y(int complex_index, int component_index,
	float x, float y, float z, float end_x, float end_y,
	float end_z, point p, float radius_sqr); // move x, z	

void calculate_stats_for_surface_z(int complex_index, int component_index,
	float x, float y, float z, float end_x, float end_y,
	float end_z, point p, float radius_sqr); // move x, y	

bool check_point(point p);

std::pair<int, int> convert_chain_index_to_atom_start_end(int chain_index);

float dist_sqr(const point& __restrict__ p1, const point& __restrict__ p2);

int get_alfa_atoms(int component_index);

int get_gridcell_index_from_point(point p);

bool is_alfa_atom(const std::string& at_name);

#endif
