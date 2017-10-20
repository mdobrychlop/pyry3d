#ifndef MUTATION_PYRY3D_H
#define MUTATION_PYRY3D_H

#include <cmath>
#include <iostream>
#include <stack>

#include "../complex.hpp"

struct anchored_history{
  char name;
  point p;
  float angle;
  char axis;
  anchored_history(char n, point p): name(n), p(p) {}
  anchored_history(char n, float an, char ax): name(n), angle(an), axis(ax) {}

};

struct grapes{
	 float radius;
	 float master_radius;
	 int grape_count;
	 std::string generation_method;

	 point master, anchor;

	 std::vector<point> list;

	 int global_index;

   std::stack< anchored_history > history;

	 bool (* detect_collisions)(int, int, const grapes&, int) = 0;


   grapes(const point& p, int g) : anchor(NAN, NAN, NAN), global_index(g){
    list.clear();
    list.push_back(p);

   }

   void copy_point(point p, int i) {
      if(i < list.size() ) {

      /* we can't copy atom index */
        list[i].x = p.x;
        list[i].y = p.y;
        list[i].z = p.z;
      }
      else
        list.push_back(p);


   }   

	 void rotate(const float angle, const char axis, bool historyb = false) {

  		float rot[27] = {1, 0, 0, 0, float(cos(angle)), float(sin(angle)), 0, //if axis is X
			float(-sin(angle)), float(cos(angle)),
				float(cos(angle)), 0, float(-sin(angle)), 0, 1, 0, 			 //if axis is Y
			float(sin(angle)), 0, float(cos(angle)),
				float(cos(angle)), float(sin(angle)), 0, 						 //if axis is Z
			float(-sin(angle)), float(cos(angle)), 0, 0, 0, 1};
  		
  	  int index; 
  		
  	  switch(axis){
  			case 'X':
  				index = 0;
  				break;
  			case 'Y':
  				index = 1;
  				break;
  			case 'Z':
  				index = 2;
  				break;
        default:
          index = -1;
  		}

  		for( int i = 0; i < list.size(); ++i ){
  			float coords[3] = {list[i].x, list[i].y, list[i].z};
  			float new_coords[3] = {0,0,0};

    		for (int j = 0; j < 3; ++j){
    					for (int k = 0; k < 3; ++k)
    						new_coords[j] += coords[k] * rot[index*9 + k*3 + j];
    			}
    			list[i].x = new_coords[0];
    			list[i].y = new_coords[1];
    			list[i].z = new_coords[2];
    		}

      if(historyb)
            history.emplace('R', angle, axis);

  		}

  		void translate(point g, bool historyb=false){
    		for(int i = 0; i < list.size(); ++i){ 
    			list[i].x += g.x;
    			list[i].y += g.y;
    			list[i].z += g.z;
    		}

        if(historyb)
            history.emplace('T', g);
  		}

      void print(int i){
        std::cout<<list[i].x<<" "<<list[i].y<<" "<<list[i].z<<"  ";
      }

      void print_history_size(){
        std::cout<<"HISTORY SIZE: "<<history.size()<<" "; 
      }
};

void exchange(int complex_index, int component_index, int func_index,
		float& angle, char& axis, std::vector<float>& l_vector, bool p,
		int& covbound_index, int atoms[]);
void rotate(int complex_index, int component_index, int func_index,
		float& angle, char& axis, std::vector<float>& l_vector, bool p,
		int& covbound_index, int atoms[]);
void rotate_cov(int complex_index, int component_index, int func_index,
    float& angle, char& axis, std::vector<float>& l_vector, bool p,
    int& covbound_index, int atoms[], int old_comp_index);
void sample(int complex_index, int chosen_indices[]);
void simul_dd(int complex_index, int component_index, int func_index,
		float& angle, char& axis, std::vector<float>& l_vector, bool p,
		int& covbound_index, int atoms[]);
void translate(int complex_index, int component_index, int func_index,
		float& angle, char& axis, std::vector<float>& l_vector, bool p,
		int& covbound_index, int atoms[]);

void mutate(int func_index, int component_index, int new_index);

float dist(const point& a, const point& b);

float get_max_angle(const point& root, const point& anchor, int i, const grapes& g);

void set_anchor_points(int complex_index, int res_index, point& p);

void start_anchored_chain();

void translate_covalent_components(int component_index, int new_index,
		std::vector<float>& l_vector, bool p);

void check_if_collide();

#endif 
