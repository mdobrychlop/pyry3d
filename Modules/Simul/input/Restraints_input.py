#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# 
# www.genesilico.pl 
#
from Modules.Trans.Interaction 		        import *

from collections					        import OrderedDict

from Modules.Simul.cpp				        import pyry3d_cpp

def flatten_logical_tree(inter):
	# there are two ararys of scores in each complex: res_scores and sym_scores
	# their indexes mean the same as indexes in the global arrays restraints and symmetries

	# now, we add third arrays, which contain logical interactions: AND and OR
	# they have indexes of interactions in restraints and symmetry arrays
	begin = pyry3d_cpp.get_restraints_size()
	end = flatten_tree(inter, begin)
	send_logical(begin, end, inter)

def flatten_tree(inter, begin):
	# postfix - firstly we want to create all interactions in restraints and
	# symmetry arrays (then we can use indexes to call them). Don't break 
	# the order, I assume the order in arrays will be the same as the order of
	# run over this tree.
	for el in inter.elements:
		if isinstance(el, SymmetryDistances):
			send_symmetry(el)
			begin+=1
		elif isinstance(el, ResidueDistanceInteraction):
			send_sd(el)
			begin+=1
		elif isinstance(el, RelationDistances):
			send_relation(el)
			begin+=1
		elif isinstance(el, SurfaceAccessInteraction):
			send_sa(el)
			begin+=1
		elif isinstance(el, PointDistanceInteraction):
			send_pd(el)
			begin+=1
		elif isinstance(el, LogicalRestraint):
			end = flatten_tree(el, begin)
			send_logical(begin, end, el)
			begin = end 
	return begin

def get_args_from_pdbranges(a, b, c, d):
	#see filtrest3D for details of classes
	loc = OrderedDict(sorted(locals().items()))
	r_list = []
	for a, i in loc.iteritems():
		r_list.append(i.chain_id)
		if isinstance(i, PDBrange):
			r_list.append(i.res1_num)
			r_list.append(i.res2_num)
		elif isinstance(i, PDBId):
			r_list.append(i.res_num)
			r_list.append(i.res_num)
		
	return r_list

def send_logical(begin, end, inter):
	pyry3d_cpp.add_logical(begin, end, isinstance(inter, AndRestraint))

def send_pd(inter):
	a = inter.restraint.reslist1[0]
	c = inter.restraint
	if not isinstance(a, PDBrange):
		pyry3d_cpp.add_point_distance(a.chain_id, a.res_num, a.res_num, \
			c.dist, c.relation, c.weight, c.point, inter.first_atom_name)
	else:
		pyry3d_cpp.add_point_distance(a.chain_id, a.res1_num, a.res2_num, \
			c.dist, c.relation, c.weight, c.point, inter.first.atom_name)

def send_relation(inter):
	args = get_args_from_pdbranges( inter.restraint.reslist1[0], \
			inter.restraint.reslist2[0], inter.restraint.reslist3[0], \
			inter.restraint.reslist4[0])
	pyry3d_cpp.add_relation_interaction(args[0], args[3], args[6], args[9], \
			args[1], args[2], args[4], args[5], args[7], args[8], args[10], \
			args[11], inter.restraint.relation, inter.first.atom_name, \
			inter.second.atom_name, inter.third.atom_name, inter.fourth.atom_name)

def send_sa(inter):
	a = inter.restraint.reslist1[0]
	c = inter.restraint
	if not isinstance(a, PDBrange):
		pyry3d_cpp.add_surface_access(a.chain_id, a.res_num, a.res_num, \
			c.dist, c.relation, c.weight, inter.first.atom_name)
	else:
		pyry3d_cpp.add_surface_access(a.chain_id, a.res1_num, a.res2_num, \
			c.dist, c.relation, c.weight, inter.first.atom_name)

def send_sd(inter):
	a = inter.restraint.reslist1[0]
	b = inter.restraint.reslist2[0] 
	c = inter.restraint
	if not isinstance(a, PDBrange) and not \
			isinstance(b, PDBrange):
		pyry3d_cpp.add_residue_distance(a.chain_id, b.chain_id, a.res_num, a.res_num, \
			b.res_num, b.res_num, c.dist, c.relation, c.weight, inter.first.atom_name, \
			inter.second.atom_name)
	elif not isinstance(a, PDBrange):
		pyry3d_cpp.add_residue_distance(a.chain_id, b.chain_id, a.res_num, a.res_num, \
			b.res1_num, b.res2_num, c.dist, c.relation, c.weight, inter.first.atom_name, \
			inter.second.atom_name)
	elif not isinstance(b, PDBrange):
		pyry3d_cpp.add_residue_distance(a.chain_id, b.chain_id, a.res1_num, a.res2_num, \
			b.res_num, b.res_num, c.dist, c.relation, c.weight, inter.first.atom_name, \
			inter.second.atom_name)
	else:
		pyry3d_cpp.add_residue_distance(a.chain_id, b.chain_id, a.res1_num, a.res2_num, \
			b.res1_num, b.res2_num, c.dist, c.relation, c.weight, inter.first.atom_name, \
			inter.second.atom_name)	

def send_symmetry(inter):
	a = inter.restraint.reslist1[0]
	b = inter.restraint.reslist2[0]
	if not isinstance(a, PDBrange) and not \
			isinstance(b, PDBrange):
		pyry3d_cpp.add_symmetry(a.chain_id, b.chain_id, a.res_num, a.res_num, \
			b.res_num, b.res_num, inter.first.atom_name, inter.second.atom_name)
	elif not isinstance(a, PDBrange):
		pyry3d_cpp.add_symmetry(a.chain_id, b.chain_id, a.res_num, a.res_num, \
			b.res1_num, b.res2_num, inter.first.atom_name, inter.second.atom_name)
	elif not isinstance(b, PDBrange):
		pyry3d_cpp.add_symmetry(a.chain_id, b.chain_id, a.res1_num, a.res2_num, \
			b.res_num, b.res_num, inter.first.atom_name, inter.second.atom_name)
	else:
		pyry3d_cpp.add_symmetry(a.chain_id, b.chain_id, a.res1_num, a.res2_num, \
			b.res1_num, b.res2_num, inter.first.atom_name, inter.second.atom_name)