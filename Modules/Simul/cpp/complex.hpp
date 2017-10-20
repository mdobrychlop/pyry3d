#ifndef PYRY3D_COMPLEX  
#define PYRY3D_COMPLEX

extern int atoms_number;
extern int chain_counter;

struct point {
	float x, y, z, atom_radii;
	int cell; //coordinates of cell in cells array -for clashes penalty computing

	point() {}

	point(float x, float y, float z): x(x), y(y), z(z), atom_radii(0) {}

	point(const std::vector<float>& f);

	point operator +=(point b){
		this->x += b.x;
		this->y += b.y;
		this->z += b.z;
		return *this;
	}

	point operator *(float b){
		point c;
		c.x = this->x * b;
		c.y = this->y * b;
		c.z = this->z * b;
		return *this;
	}

	point operator /(float b){
		point c;
		c.x = this->x / b;
		c.y = this->y / b;
		c.z = this->z / b;
		return *this;
	}

	bool operator== (const point& other) const {
		// if coordinates are the same, we can assume that it is the point we are considering
		// don't overuse this operator
		return x == other.x && y == other.y && z == other.z;
	}

	point reversed_point(){
		/* returns point symetrical by 0,0,0 point */
		point temp;
		temp.x = -x;
		temp.y = -y;
		temp.z = -z;
		temp.atom_radii = atom_radii;
		return temp;
	}
};


std::ostream& operator<<(std::ostream& os, const point& dt);

struct grid {
/// radius is half of gridCell size
///it is with overlap
	float sb_xmin;
	float sb_xmax;
	float sb_ymin;
	float sb_ymax;
	float sb_zmin;
	float sb_zmax;
	float xmin;
	float xmax;
	float ymin;
	float ymax;
	float zmin;
	float zmax;
	float radius; //cały bok siatki
	std::vector<float> grid_cells_density;
	std::vector<bool> grid_cells_ismap;
	int mapcells;//gridcell with map inside
	float density_sum;
};

struct component {
	float centre_of_mass[3];
	
	/* 3 axes + 'L' */
	
	float rot_ranges_sum[4][2];	
	float rot_history_sum[4];	

	/* 3 axes */
	float trans_history_sum[3];
	float trans_ranges_sum[3][2];
	
	bool changed = false;
	bool moves_limited = false;
	
	//scores
	int alfa_atoms = 0;
	bool is_outmap = false;
	bool is_outbox = false;
	int outbox_atoms = 0;
	int outbox_cells = 0;
	int outmap_cells = 0;
	int mapcells_nr = 0;
	float taken_density = 0.0;
	std::vector<int> taken_mapcells;
	uint32_t counted_taken_mapcells = 0;
	std::vector<float> sphere_dist_penalties;
};

struct complex {
	
	std::vector<uint16_t> Grid_taken_mapcells;

	//coordinates
	std::vector<point> coords;

	//represents clashes between two indexed components
	std::vector<std::vector<int> > clashes_nr;

	//represents clashes between sphere
	std::vector<std::vector<int> > sphere_penalties;
	
	//current score
	float score = 0;

	//last remembered iteration
	int iteration_index;
	
	//temperature
	float reptemp, temp;

	//these are for counting symmetry and restraints scores
	std::vector<float> sym_scores, res_scores;
	std::vector<bool> already_visited_by_logical;
	std::vector<bool> rearranged;

	float symmetry, restraints;
	
	//number of alfa atoms in complex
	int alfa_atoms;

	float all_clashes;
	int all_clashes_dist;
	int all_filled_cells;
	int all_outbox_atoms;
	float all_restraints;
	float clashes;
	float density;
	float freespace;
	bool is_complex_outbox;
	bool is_complex_outmap;
	int mapcells_nr;
	float outbox;
	int pairs_atoms;
	float taken_densities;
	int taken_mapcells;
	float rg;
	float chi2;
	
	std::vector<component> components;
	
	complex();

	complex(std::vector<point> p): coords(p), score(0.0f), temp(0.0f) {}

	complex(std::vector<point> p, std::vector<std::vector<int> > vec_clashes_nr,
			std::vector<std::vector<int> > vec_all_clashes_dist):
			coords(p), clashes_nr(vec_clashes_nr), sphere_penalties(vec_all_clashes_dist),
			score(0.0f), temp(0.0f) {}


};

struct simul_params {
	bool required_clashes_penalty = false;
	bool required_clashes_penalty_allatoms = false;
	float outbox[2];
	float freespace[2];
	float clashes[2];
	float restraints[2];
	float density[2];
	float symmetry[2];
	float chi2[2];
	float rg[2];
	float rg_val;
	std::string curve;
	std::string crysol_path;
	std::string name_prefix;
	std::string representation;
};

struct by_simulation_score_descending { 
    const bool operator()(complex const &a, 
		complex const &b) const ;
};

#endif

