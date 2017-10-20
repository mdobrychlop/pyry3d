#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# 
# www.genesilico.pl 
#
import os

from Modules.Error.Errors                   import SimulError
from math							        import exp, pow
from random							        import random
from Bio.PDB.Model                          import Model
from Bio.PDB                                import PDBIO
from Bio.PDB.Structure						import Structure
from Modules.Constans.Logfile				import outfolder
from Modules.Constans.Constans				import ATOM_RADII
from Modules.Error.Errors					import InputError
from Modules.Simul.cpp						import pyry3d_cpp
from Modules.Simul.input.Restraints_input	import *
from Modules.Simul.Test 					import *
from Modules.Input.Input					import Input	
from Modules.Simul.Mutator					import Factory
from Modules.Input.DataDownloader           import Sequences, Structures
from Modules.Trans.Trans                    import Trans
from Modules.Simul.Complex                  import PyRyComplex
from Modules.Simul.PlotComplex                 	import PlotComplex
from Modules.Simul.Send_first_complex      import send_first_grid, send_first_complex, send_first_clashes, send_coords
from Modules.Config.Config                  import Config
from Modules.Constans.Constans		        import MOLWEIGHTS
from copy                           import deepcopy
from External_Applications.filtrest3D.PDBId import PDBId
from External_Applications.filtrest3D.PDBrange import PDBrange

from sys 							import exit

def temperature (step, aTemp):
	"""
	returns new temperature value
	"""
        alpha = 0.999
        return aTemp * pow(alpha, step)

def transfer_to_pyrycomplex(i, sim):
	"""
	returns pyry Complex from cpp module (only coords)
	"""
	res_complex = deepcopy(sim.mPool[0])
	coords = pyry3d_cpp.get_coords(i)
	i = 0;
	for component in res_complex.components:
		for at in component.pyrystruct.struct.get_atoms():
			at.coord[0] = coords[i][0];
			at.coord[1] = coords[i][1];
			at.coord[2] = coords[i][2];
			i += 1
	
	return res_complex

def cpp_anneal_complex(index, iteration, max_temperature, sim):
	# New index is index of new created complex in complexes array;
	# Note that new index is out of pool 
	new_index = pyry3d_cpp.perform_mutation(index, iteration)
	cpp_calculate_rating(new_index, sim)
	current_temp = temperature(iteration, max_temperature)
	pyry3d_cpp.set_temp(new_index, current_temp)  
	if ( cpp_metropolis_compare(max_temperature, new_index, index, sim, current_temp) ):
		# Accept the change
		pyry3d_cpp.exchange(new_index, index)
	# Remove last complex (we don't need it)
	pyry3d_cpp.pop()
	
def cpp_calculate_rating(index, sim):

	"""
	equivalent to calculateRating from Simul
	"""

	pyry3d_cpp.update_restraints_score(index)
	pyry3d_cpp.calculate_restraints_score(index)
	pyry3d_cpp.calculate_symmetry_score(index)
	
	if pyry3d_cpp.number_of_components() == 1:
		pyry3d_cpp.assign_new_taken_indexes(index)
		pyry3d_cpp.calculate_score_for_one_component_complex(index)
		
	else:
		pyry3d_cpp.set_false_is_complex_outbox(index)
		pyry3d_cpp.set_false_is_complex_outmap(index)
		pyry3d_cpp.assign_new_taken_indexes(index)
		pyry3d_cpp.calculate_outbox_mapfill(index)
		pyry3d_cpp.calculate_simulation_score_for_complex(index)	
		
def cpp_choose_result_save_procedure(iteration, sim):	
	
	if sim.simul_params.save_res_mode == "eachbetter":
		if pyry3d_cpp.get_simulation_score(0) > sim.lBestSimulation_score:
			cpp_save_simulation_results(iteration, [0], sim)
			sim.lBestSimulation_score = pyry3d_cpp.get_simulation_score(0)
		
	if sim.simul_params.movehistory and pyry3d_cpp.get_simulation_score(0) \
			> pyry3d_cpp.get_best_complex_score():
		# copies complex at 0 index to best_score_complex variable
		pyry3d_cpp.deepcopy_complex_to_best_score_complex(0)
		
	if sim.simul_params.save_res_mode == "niter":
		if (iteration % sim.mTrajStep == 0):
			if (sim.mSimulationMethod == sim.SimulationMethods.ReplicaExchange):
				# saves whole pool
				cpp_save_simulation_results(iteration, range(0, pyry3d_cpp.pool_size()), sim)
			outname = outfolder.outdirname.split("/")[-1]
			name = str(outfolder.outdirname)+'/'+str(outname)+"_"
			pyry3d_cpp.save_best_niter(iteration, name)
			# copies complex at 0 index to best_score_complex_niter variable
			pyry3d_cpp.deepcopy_complex_to_best_score_niter_complex(0)
	    
		else:
			if (pyry3d_cpp.get_simulation_score(0) > pyry3d_cpp.get_niter_score()):
				pyry3d_cpp.deepcopy_complex_to_best_score_niter_complex(0)
	
	if sim.simul_params.save_res_mode == "outsteps":
	    cpp_collect_complexes_to_save(iteration, 0, sim)
	    
	#here we increment step
	pyry3d_cpp.inc_step()

def cpp_collect_complexes_to_save(iteration, complex_id, sim):
	if sim.simul_params.out_steps[0] <= iteration <= sim.simul_params.out_steps[1]:
		if pyry3d_cpp.number_of_to_save() < sim.simul_params.struct_nr:
			# if we remember less than struct_nr complexes, we just add another one
			pyry3d_cpp.add_to_save(complex_id)
		else:
			# else we need to check if there is complex that we can replace
			pyry3d_cpp.try_to_push_in_to_save(complex_id)
	
def cpp_init(sim):

	# 90% of code here is copying data from python to c++

	disorder_config_check(sim)

	#we can expect here different type - can't be used this way in C++, so
	#we need an ugly solution
	if isinstance(sim.simul_params.shapedesc.__class__, bool):
		shapedesc = sim.simul_params.shapedesc
	else:
		shapedesc = True

	# The following sums get number of atoms (don't be afraid)
	pyry3d_cpp.init_cpp(sum(sum(1 for x in comp.pyrystruct.struct.get_atoms()) \
			for comp in sim.mPool[0].components), sim.mPool[0].clashes_penalty,\
			sim.mPool[0].restraints_penalty, sim.mPool[0].freespace_penalty, \
			sim.mPool[0].outbox_penalty, sim.mPool[0].density_penalty,\
			sim.mPool[0].symmetry_penalty, sim.mPool[0].chi2_penalty,\
			sim.mPool[0].rg_penalty, sim.mIterationCount, \
			sim.scaling == "on", len(sim.mPool[0].components), MOLWEIGHTS, \
			shapedesc, str(outfolder.outdirname)+'/')

	pyry3d_cpp.get_transformation_frequencies([sim.rotation_freq, \
			sim.translation_freq, sim.exchange_freq, sim.simul_dd, \
			sim.translation_all_freq, sim.rotation_all_freq, \
			sim.rotation_cov_freq, sim.rotation_whole_freq, sim.exchangeandsample])

	pyry3d_cpp.send_arrays( sim.mPScalRangeBounds ,sim.mPScalValues,\
			sim.mPool[0].movable, sim.mPool[0].free, \
			[a.coord for a in sim.mDensityMap.surfacepoints])

	temp = sim.simul_params
	outname = outfolder.outdirname.split("/")[-1]

	curve = "" if temp.curve == None else temp.curve
	crysol_path = "" if temp.crysol_path == None else temp.crysol_path

	pyry3d_cpp.get_simul_params(temp.clashes_penalty[0],\
			temp.clashes_penalty[1], temp.restraints_penalty[0],\
			temp.restraints_penalty[1], temp.freespace_penalty[0],\
			temp.freespace_penalty[1], temp.outbox_penalty[0],\
			temp.outbox_penalty[1], temp.density_penalty[0],\
			temp.density_penalty[1], temp.symmetry_penalty[0], \
			temp.symmetry_penalty[1], temp.chi2_penalty[0], \
			temp.chi2_penalty[1], temp.rg_penalty[0], \
			temp.rg_penalty[1], curve, crysol_path, \
			str(outfolder.outdirname)+'/'+str(outname)+"_", \
			temp.rg_val, str(temp.representation), \
			temp.required_clashes_penalty_allatoms, temp.required_clashes_penalty)
	

	# file has to be send before we continue with components
	send_coords(sim.mPool[0],0)
	
	for c in sim.mPool[0].components:
		#here we copy components data (disorders, covalent bonds, etc.)

		pyry3d_cpp.create_component(c.moves.allowed_transform, c.mass_centre,
				c.moves.trans_ranges_sum, c.moves.trans_ranges, c.moves.limited,
				c.moves.rot_axis, c.moves.rot_ranges, c.moves.rot_ranges_sum,
				c.pyrystruct.moltype, len(c.disorders) > 0,
				len(sim.mPool[0].components))
		for comp_index in c.covalent_bonds:
			pyry3d_cpp.add_covalent_bond(comp_index.chains_indexes, \
					comp_index.atom1[0], comp_index.atom2[0]) 
		for disorder in c.disorders:
			pyry3d_cpp.add_disorder(disorder.start_pos, disorder.stop_pos, \
					disorder.fragment_type, disorder.radius, 
					disorder.max_sphere_radius)
	
	#copying restraints

	for inter in sim.mInteractions.symmetry_interactions:
		send_symmetry(inter)
	
	for inter in sim.mInteractions.sd_interactions:
		send_sd(inter)

	for inter in sim.mInteractions.sa_interactions:
		send_sa(inter)
	
	for inter in sim.mInteractions.pd_interactions:
		send_pd(inter)

	for inter in sim.mInteractions.relation_interactions:
		send_relation(inter)

	# those have to be added after all other interactions - order matters!

	# translates recursion in restraints (And, Or restraints) to iteration
	for inter in sim.mInteractions.logical_interactions:
		flatten_logical_tree(inter)

	# creates new array of component indexes which can be exchanged
	pyry3d_cpp.fill_can_be_exchanged()

	#sends grid used in density penalty calculation
	send_first_grid(sim.mDensityMap)

	sphere_radii = []
	if sim.simul_params.representation == "sphere":
		for comp in sim.mPool[0].components:
			for atom in comp.pyrystruct.struct.get_atoms():
				sphere_radii.append(atom.vdw)

	pyry3d_cpp.send_atom_radii(ATOM_RADII, sphere_radii)
	send_first_complex(sim.mPool[0], 0) # all data for first complex and its components
	send_first_clashes(sim.mPool[0], 0) # see above	

	pyry3d_cpp.calculate_mass_of_complex_and_componets(0);
			

def cpp_metropolis_compare(mMaxTemperature, complex_a, complex_b, sim, current_temp = 0 ):
	if current_temp == 0:
		current_temp = mMaxTemperature
		
	diff = pyry3d_cpp.get_simulation_score(complex_a) - pyry3d_cpp.get_simulation_score(complex_b)
	if current_temp < 0.0001:
		current_temp = 0.0001

	acceptance_probability = 0 if ( (-diff/current_temp) > 700.0) \
		else 1.0 / (1.0 + exp( (-diff)/current_temp ) )
	
	rand = random.random()

	if ( (diff > 0) or (rand < acceptance_probability) ):    #Accept the change
		print "CHANGE ACCEPTED!!!"
		sim.x += 1
		return True
		
	return False

def cpp_save_last_simulated_complex(iteration, sim, sim_methods):
	
	if sim.simul_params.save_res_mode == "outsteps":
		pyry3d_cpp.try_to_push_in_to_save(-1)
		#with -1 as an argument it only sorts complexes
		outname = outfolder.outdirname.split("/")[-1]
		name = str(outfolder.outdirname)+'/'+str(outname)+"_"

		logfile_data = pyry3d_cpp.save_complexes_to_save(sim.mReplicaTemperatures != [], name)
		
		for i in range(0, len(logfile_data), 7):
			
			replica = createReplicaValuesOutstep(i, logfile_data, sim)
			
			if (replica.clashes == -1):
				logfile.write_file("Rejected\t" + str(replica.simulation_score) +\
				"\titer_nr\t" + str(logfile_data[i] + 1) + "\n")
			else:
				logfile.write_file(\
				 "outsteps cx score for iteration "+str(logfile_data[i] + 1)+" is\t" + str(round(replica.simulation_score, 3))\
				+"\tcomponents: restraints: "+str(round(replica.restraints, 3) * sim.restraints_penalty)+" "+\
				"collisions: "+str(round(replica.clashes, 3) * sim.collision_penalty)+" "+\
				"map filling: "+str(round(replica.freespace, 3) * sim.freespace_penalty)+" "+\
				"outbox atoms: "+str(round(replica.outbox, 3) * sim.outbox_penalty)+" "+\
				"density filling: "+str(round(replica.density, 3) * sim.density_penalty)+" "+\
				"\tACTUAL WEIGHTS are respectively: "+str(sim.restraints_penalty)+", "+\
				str(sim.collision_penalty)+", "+str(sim.freespace_penalty)+", "+\
				str(sim.outbox_penalty)+", "+str(sim.density_penalty)+"\n")

				if sim.plot:
					plot_cx = PlotComplex(replica, logfile_data[i])
					sim.saved_complexes.append(plot_cx)
				
				print "should be saved " + str(logfile_data[i])
	
			
	if sim.simul_params.save_res_mode == "niter":
		
		if (sim.mSimulationMethod == sim_methods.ReplicaExchange):
			cpp_save_simulation_results(iteration, \
					range(0, pyry3d_cpp.pool_size()), sim)
		else: 
			cpp_save_simulation_results(iteration, [0], sim)
	
def cpp_save_simulation_results(iteration, complexes, sim):
		
	for c in complexes:
		if not pyry3d_cpp.is_in_scores(c):
			pyry3d_cpp.append_score(c)
			outname = outfolder.outdirname.split("/")[-1]
			name = str( outfolder.outdirname )+'/'+str( outname )+"_"+\
					str( round(pyry3d_cpp.get_simulation_score(c), 4 ) )+'_'+\
					str( pyry3d_cpp.get_step() )+"_"
			if sim.mReplicaTemperatures:
				name = name+str(pyry3d_cpp.get_reptemp(c) )+'.pdb'
			else:
				name = name+str(pyry3d_cpp.get_temp(c) )+'.pdb'
			pyry3d_cpp.write_to_file(c, name, pyry3d_cpp.get_atoms_number() ) 
			
			logfile.write_file(\
		    "cx score for iteration " + str(iteration + 1) + " is\t"\
		    + str(round(pyry3d_cpp.get_simulation_score(c), 3))\
		    +"\tcomponents: restraints: "+str(round(pyry3d_cpp.get_simulation_restraints(c), 3) * sim.restraints_penalty)+" "+\
		    "collisions: "+str(round(pyry3d_cpp.get_simulation_clashes(c), 3) * sim.collision_penalty)+" "+\
		    "map filling: "+str(round(pyry3d_cpp.get_simulation_freespace(c), 3) * sim.freespace_penalty)+" "+\
		    "outbox atoms: "+str(round(pyry3d_cpp.get_simulation_outbox(c), 3) * sim.outbox_penalty)+" "+\
		    "density filling: "+str(round(pyry3d_cpp.get_simulation_density(c), 3) * sim.density_penalty)+" "+\
		    "\tACTUAL WEIGHTS are respectively: "+str(sim.restraints_penalty)+", "+\
		    str(sim.collision_penalty)+", "+str(sim.freespace_penalty)+", " +\
		    str(sim.outbox_penalty)+", "+str(sim.density_penalty)+"\n")

			if sim.plot:
				replica = createReplicaValues(c, sim)
				plot_cx = PlotComplex(replica, iteration)
				sim.saved_complexes.append(plot_cx)
				
		else:
			logfile.write_file("Rejected\t" + str(pyry3d_cpp.get_simulation_score(c)) +\
			"\titer_nr\t" + str(iteration + 1) + "\n")
	
	sim.saved_iterations.append(iteration)


def createReplicaValues(complex_index, sim):

	"""
	returns PyRyComplex with data for plot
	"""
	replica = PyRyComplex()
	replica.simulation_score = pyry3d_cpp.get_simulation_score(complex_index)
	replica.clashes = pyry3d_cpp.get_simulation_clashes(complex_index)
	replica.clashes_penalty = sim.mPool[0].clashes_penalty
	replica.density = pyry3d_cpp.get_simulation_density(complex_index)
	replica.density_penalty = sim.mPool[0].density_penalty
	replica.freespace = pyry3d_cpp.get_simulation_freespace(complex_index) 
	replica.freespace_penalty = sim.mPool[0].freespace_penalty
	replica.outbox = pyry3d_cpp.get_simulation_outbox(complex_index) 
	replica.outbox_penalty = sim.mPool[0].outbox_penalty
	replica.restraints = pyry3d_cpp.get_simulation_restraints(complex_index) 
	replica.restraints_penalty = sim.mPool[0].restraints_penalty
	
	return replica

def createReplicaValuesOutstep(index, logfile_data, sim):

	"""
	returns PyRyComplex with data for save_res_mode == "outsteps"
	"""
	replica = PyRyComplex()	
	replica.simulation_score = logfile_data[index + 1]
	replica.clashes = logfile_data[index + 3]
	replica.clashes_penalty = sim.mPool[0].clashes_penalty
	replica.density = logfile_data[index + 6]
	replica.density_penalty = sim.mPool[0].density_penalty
	replica.freespace = logfile_data[index + 4]
	replica.freespace_penalty = sim.mPool[0].freespace_penalty
	replica.outbox = logfile_data[index + 5] 
	replica.outbox_penalty = sim.mPool[0].outbox_penalty
	replica.restraints = logfile_data[index + 2]
	replica.restraints_penalty = sim.mPool[0].restraints_penalty

	return replica


def disorder_config_check(sim):
	factory = Factory(sim.configuration, sim)
	if(len(sim.mPool[0].with_disorders) == 0 and
			"SimulateDisorder" in factory.configuration.mutation_frequencies and
			factory.configuration.mutation_frequencies["SimulateDisorder"] > 0):
		raise InputError("You cannot apply Simulate Disorder mutation for " +
			"components with no flexible/disordered fragments")