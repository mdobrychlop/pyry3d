#ifndef RESTRAINTS_PYRY3D 
#define RESTRAINTS_PYRY3D 

#include <vector>

#include "../complex.hpp"

class ResidueDistance;
class SurfaceAccess;

class Restraint{
	public:

		/* types:
			's' - symmetry
			'd' - distance
			'r' - relation with distance
			'q' - relation with two pairs of residues
		*/

		char type;
		int first;
		int second;

		int first_range[2];
		int second_range[2];

		std::string first_atom_name;
		std::string second_atom_name;

		void score(int complex_index, int index, int first, int second,
			int third, int fourth);
		virtual void get_score(int complex_index, int index, int first, int second,
			int third, int fourth) = 0;
		bool operator==(const Restraint &q);
		virtual ~Restraint() {}
		Restraint();
		Restraint(int a1, int a2, int b1, int b2, char t, const std::string& f_a_n, 
				const std::string& s_a_n);
};

class Logical {
	public:
		int id;		/* index in array */
		int begin, end;
		bool max; 
		/* if yes then and if not then or  - i prefer max/min notation*/
		void score(int complex_index);
		Logical(int b, int e, bool and_);
};

class RelationRestraint : public Restraint{
	/* relation 1 - "<" 2 - "<=" 3 - ">=" 4 - ">" */
	public:
		short relation;
		float distance;	/* ex. "<=" (4.098) */
		float weight;
		void print();
		RelationRestraint(int a1, int a2, int b1, int b2, float dist, 
				const std::string& rel, float weig, char t, const std::string& f_a_n,
				const std::string& s_a_n);
		virtual ~RelationRestraint() {}
		virtual void get_score(int complex_index, int index, int first, int second) = 0;
		void get_score(int complex_index, int index, int first, int second,
			int third, int fourth);
};

int calculate_diff(int complex_index, int i, int j);
int simplify_relation(int rel);

class PointDistance : public RelationRestraint{
	public:
		point p;

		void get_score(int complex_index, int index, int first, int second);

		PointDistance(const std::string& s,	int a1, int a2,	float dist, 
				const std::string& rel, float weig, 
				const std::vector<float>& point_v, const std::string& f_a_n);
};

class ResidueDistance : public RelationRestraint{
	public:

		void get_score(int complex_index, int index, int first, int second);

		ResidueDistance(const std::string& s, const std::string& f,  
				int a1, int a2, int b1, int b2,
				float dist, const std::string& rel, float weig, 
				const std::string& f_a_n, const std::string& s_a_n);

};

class RelationDistance : public Restraint{
	private:
		void get_pair_score(int complex_index, int first, int second,
				int f_range[], int s_range[], float& min_diff, float& max_diff,
				const std::string& f_name, const std::string& s_name);

	public:
		short relation;

		int third;
		int fourth;

		int third_range[2];
		int fourth_range[2];

		std::string third_atom_name;
		std::string fourth_atom_name;

		void get_score(int complex_index, int index, int first, int second, int third,
				int fourth);

		RelationDistance(const std::string& f, const std::string& s,
				const std::string& t, const std::string& fo,  
				int a1, int a2, int b1, int b2, int c1, int c2, int d1, int d2,
				const std::string& rel, const std::string& f_a_n, 
				const std::string& s_a_n, const std::string& t_a_n,
				const std::string& fo_a_n);

};

class SurfaceAccess : public RelationRestraint{
	public:

		void get_score(int complex_index, int index, int first, int second);

		SurfaceAccess(const std::string& s, int a1, int a2,
				float dist, const std::string& rel,	float weig,
				const std::string& f_a_n);
};



void calculate_restraints_score(int index);
void update_restraints_score(int index);

#endif 