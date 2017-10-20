#ifndef PYRY3D_SEND_DATA_HPP  
#define PYRY3D_SEND_DATA_HPP  

extern std::vector<complex> complexes;

void calculate_first_complex_taken_mapcells_stats(int complex_index, int i);

//~ void send_first_pdb(int complex_index, std::string send_pdb);

void set_changed_component(int complex_index, const std::vector<int>& changes);
void set_clashes_nr(int complex_index, std::vector<int> clashes);
void set_complex_data(int complex_index, int alfa_atoms,
	float all_clashes, int all_filled_cells, int all_outbox_atoms,
	float all_restraints, float clashes, float density,
	bool is_complex_outbox, bool is_complex_outmap, int mapcells_nr,
	int pairs_atoms, float restraints, float simulation_score,
	float symmetry, float taken_densities, int taken_mapcells);
	
void set_components_data(int complex_index, std::vector<bool> is_outbox,
	std::vector<bool> is_outmap, std::vector<int> mapcells_nr,
	std::vector<int> outbox_atoms);

#endif
