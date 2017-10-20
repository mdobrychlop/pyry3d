#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# 
# www.genesilico.pl 
#

#creates ranked 3D models of macromoleular complexes 
#based on experimental restraints and a whole complex shape.

__author__ = "Joanna M. Kasprzak"
__copyright__ = "Copyright 2010, The PyRy3D Project"
__credits__ = ["Janusz Bujnicki"]
__license__ = "GPL"
__version__ = "0.1.0"
__maintainer__ = "Joanna Kasprzak"
__email__ = "jkasp@amu.edu.pl"
__status__ = "Prototype"

import time, warnings
#NumPy:
from numpy                          import array
#standard libraries:
from random                         import uniform, seed
from shutil                         import rmtree
import os
from copy                           import deepcopy
#imports from PyRy3D:
from Modules.Input.Input            import Input
from Modules.Config.Config          import Config
from Modules.Trans.Trans            import Trans
from Modules.Simul.Simul            import Simul
from Modules.Simul.Complex          import PyRyComplex
from Modules.Constans.Logfile       import logfile, trafl, movehist, outfolder
from Modules.Error.Errors           import InputError
#
##plot draw libs
#~ from pylab import *
#~ import matplotlib.pyplot as plt
#~ import matplotlib.axis as axi
#~ import matplotlib.ticker as ticker


def arrange_components_in_map_center(components, density_map, direct="ahead"):
    """
        Puts all components inside simulbox;
        After rearangment all centers of mass (for all complex components) are localized inside a cuboid center which represents density map space
        This method is called for statistical purposes only;
        It allows to calculate number of grid cells necesarry to describe each component.
    Parameters:
    ----------
        components  : list of all complex components
        density_map : density map object
        direct      : "ahead" for first translation, "reverse" to arrange components back to start position
    """
    for component in components:
        component.pyrystruct.calculate_centre_of_mass()
        vector = []
        if direct == "reverse":
            vec = component.trans_history
            component.translate(-array(vec[1]))
        elif direct == "ahead":
            vector = array(density_map.mapcuboid.center) - array(component.pyrystruct.center_of_mass)
            component.translate(vector)        
                

def arrange_components_inside_simulbox(fcomplex, density_map, trans, con):
    """
        puts all components inside the simulation box;
        after rearangment all centers of mass are randomly localized inside a density map
    Parameters:
    ----------
        components  : list of all complex components
        density_map : density map object
        trans and con are PyRy3D objects from Trans.py and Config.py modules
    """
    
    print "PyRy3D arranged complex components randomly inside the SIMULATION BOX"
    
    print "The PyRy3D score now is equal to\t:"
    
    if con.start_orient == True:
        #pass
        for component in fcomplex.components:
            if len(component.disorders) != 0:
                for el in component.disorders:
                    if el.fragment_type == "simulated_volume":
                        print "Component with unknown structure will be positioned randomly inside simulation box"
                        position_component(component)
                        component.is_positioned = True
    else:
        for component in fcomplex.components:
            if component.moves.state == "fixed" or component.moves.limited == True:
                #do not move components randomly if the user provided any restraints on their moves eg fixed or move_state are provided
                pass
            elif len(component.covalent_bonds) != 0:
                #position this anchor component randomly in the space
                if component.is_positioned == False:
                    vec = position_component(component)
                    #position linked components by translating them by the exact the same vector
                    for covbond in component.covalent_bonds:
                        for comp_index in covbond.chains_indexes:
                            comp = fcomplex.components[comp_index]
                            print "components for covbond", comp.pyrystruct.chain, vec
                            comp.translate(vec)
                            comp.is_positioned = True
            else:
                if component.is_positioned == False:
                    position_component(component)
                    component.is_positioned = True
    
    fcomplex.calculate_simulation_score(trans.interactions, density_map, con, iteration="first_complex")
    fcomplex.update_restraints_score(trans.interactions)
    fcomplex.save_pdb('inside_map')
    fcomplex.calculate_centre_of_complex()


    
def position_component(component):
    """
    Method translates a component from start orientation (as provided in .pdb file)
    into random orientation inside a density map
    """
    component.pyrystruct.calculate_centre_of_mass()
    random_vector = calculate_translation_vector(density_map)
    vector = array(random_vector) - array(component.pyrystruct.center_of_mass)
    component.translate(vector)
    return vector
            
def calculate_translation_vector(density_map):
    """
        calculates how to translate particular component in order
        to locate it inside the mapcuboid new position is chosen randomly
    Parameters:
    ----------
        density_map :   density map Complex_map object
    Returns:
    --------
        random_vector   : randomly chosen [x,y,z] coordinates inside mapcuboid
    """
    random_x = uniform(density_map.mapcuboid.xmin, density_map.mapcuboid.xmax)
    random_y = uniform(density_map.mapcuboid.ymin, density_map.mapcuboid.ymax)
    random_z = uniform(density_map.mapcuboid.zmin, density_map.mapcuboid.zmax)
    random_vector = [random_x, random_y, random_z]
    return random_vector
        
def calculate_simulbox_statistics(fcomplex, density_map,inp, trans, con):
    """
    if program works with complex shape a function puts components into its center in order to calculate
    scoring function statistics; in next step components are randomly oriented in simulation box
    """
    
#@TODO:this is not creation of first complex but seting map and calculation of scores!! 
    fcomplex.save_pdb('original_structure')    
      
    if inp.mapfile or inp.saxsfile:
    #calculation of statistics for simbox cells occupation is neccesary only when shape is defined!
        arrange_components_in_map_center(fcomplex.components, density_map)
        fcomplex.save_pdb('incenter')
        fcomplex.get_simboxcells(density_map)
        fcomplex.save_pdb('simboxcells')
        arrange_components_in_map_center(fcomplex.components, density_map, "reverse")
        fcomplex.save_pdb('reversed')
    
    print "input complex score:"
    fcomplex.calculate_simulation_score(trans.interactions, density_map, con, iteration="first_complex")
    
    
def check_corretness_of_covalent_bonds(con):
    """
    check whether all chains defined in covalent bonds exist in provided structures
    here all components chains names must exist in pdb files provided
    if some components do not have covalent bonds defined they are treated as independent and are mutated alone.
    """
    for chain in con.linked_components.keys():
        if chain not in con.components_indexes.keys():
            raise InputError("chain %s does not occur in structures you provided. please correct covalent bonds list in configuration file"%(chain))
        for el in con.linked_components[chain]:
            for e in el.chains:
                if e not in con.components_indexes.keys():
                    raise InputError("!!chain %s does not occur in structures you provided. please correct covalent bonds list in configuration file"%(e))

def create_first_complex(components, con, traflfile):
    """
        Sets the initial structure components start_anntemp :
        starting simulation temperature
    Parameters:
        components  : structures in the complex
        con         : PyRy3D object from Config.py
        traflfile   : output file object storing simulation trajectory
    """
#@TODO: clean this function!!!

#----- call first complex, assign components and start temperature ---------
    pc = PyRyComplex (components, None)
    
    
    pc.set_annealing_temp(con.anntemp)
    pc.set_penalties_weights(con.clashes_penalty[0],\
                             con.restraints_penalty[0],\
                             con.freespace_penalty[0],\
                             con.outbox_penalty[0],\
                             con.density_penalty[0],\
                             con.symmetry_penalty[0],\
                             con.chi2_penalty[0],\
                             con.rg_penalty[0])
    
# -------------------- add first complex to outfiles ----------------------    
    logfile.write_message("cx score for iteration 0 is\t"+str(round(pc.simulation_score,3))+"\tcomponents: "+\
                         "restraints: "+str(round(pc.restraints,3) * con.restraints_penalty[0])+" "+\
                         "collisions: "+str(round(pc.clashes,3)    * con.clashes_penalty[0])+" "+\
                         "map filling: "+str(round(pc.freespace,3)* con.freespace_penalty[0])+" "+\
                         "outbox atoms: "+str(round(pc.outbox,3)     * con.outbox_penalty[0])+" "+\
                         "density filling: "+str(round(pc.density,3)  * con.density_penalty[0])+\
                         "symmetry: "+str(round(pc.symmetry,3)  * con.symmetry_penalty[0])+\
                         "chi2: "+str(round(pc.chi2,3)  * con.chi2_penalty[0])+\
                         "RGE: "+str(round(pc.rg,3)  * con.rg_penalty[0])+\
                         "\tACTUAL WEIGHTS are respectively: "+str(con.restraints_penalty[0])+", "+\
                         str(con.clashes_penalty[0])+", "+str(con.freespace_penalty[0])+", "+\
                         str(con.outbox_penalty[0])+", "+str(con.density_penalty[0])+", "+str(con.symmetry_penalty)+", "+\
                         str(con.chi2_penalty)+", "+str(con.rg_penalty))

# -------------------------------------------------------------------------

#---- set complex properties before first simulation run -------
    
    for component in pc.components:
        component.clear_pyryatoms()
        component.simulate_disorder(component.pyrystruct.struct)
        component.add_covalent_bonds(con, pc)
       
#@TODO:##these might be combined into one function in Complex class
        pc.add_simul_change(component)
        
        pc.add_movable_component(component)
       
        pc.get_alfa_atoms(component)
       
        pc.add_complex_atoms(len(list(component.pyrystruct.struct.get_atoms())))
        
        pc.add_complex_alfa_atoms(len(component.alfa_atoms))
    pc.calculate_complex_volume ()
    pc.generate_comp_pairs()
    pc.clean_vollist()
    pc.diffscores = []
    return pc

def draw_scores_plot(plot_name, sim):
    """
    method enables to draw energy plots for all complexes generated during simulation process;
    As an output it returns two plots:
    scoreelems.plot : with all scoring function elements values
    simscores.plot  : just with general (total) score values
    Parametrs:
       plot_name : name of the output plot
       sim       : Simul.py object
    """
    saved_complexes = sim.get_plot_values()
    
    scores, clashes, restraints, outbox, mapfill, densities, symmetries, steps = [],[],[],[],[],[],[], []
    for cx in saved_complexes:
        scores.append(cx.sim_score)
        clashes.append(cx.clashes)
        restraints.append(cx.restraints)
        outbox.append(cx.outbox)
        mapfill.append(cx.mapfill)
        densities.append(cx.densities)
        symmetries.append(cx.symmetry)
        steps.append(cx.step_nr)
        
        
    print "Scores", scores
    #Rysowanie wtkresu
    fig = plt.figure()
    xlabel("Steps", fontsize='x-large')
    ylabel("Scores", fontsize='x-large')
    #plot(steps, scores, color ='black', marker = '.', ms=8)
    plot(steps, clashes, color ='green', marker = '.', ms=8)
    plot(steps, restraints, color ='blue', marker = '.', ms=8)
    plot(steps, outbox, color ='red', marker = '.', ms=8)
    plot(steps, mapfill, color ='cyan', marker = '.', ms=8)
    plot(steps, densities, color ='pink', marker = '.', ms=8)
    plot(steps, symmetries, color ='yellow', marker = '.', ms=8)
    plt.legend(("clashes","restraints","outbox", "density map filling","density filling", "symmetry"))
    #show()
    fig.savefig(str(plot_name)+"_scoreelems.plot", dpi = 500, format = 'png', transparent = False)
    
    fig2 = plt.figure()
    xlabel("Steps", fontsize='x-large')
    ylabel("Scores", fontsize='x-large')
    plot(steps, scores, color ='black', marker = '.', ms=8)
    fig2.savefig(str(plot_name)+"_simscores.plot", dpi = 500, format = 'png', transparent = False)
   
def get_input_data(inp, con, trans):
    """
        calls Input Module methods to get data delivered by the user i.e:
            restraints_file  :   file with distance restraints in Filtrest3D format
            structures_files :   tar archive with all pdb structures;
                                 ATTENTION! one file can store one chain only
            config_file      :   file with simulation parameters defined by the user
            sequence_file    :   plain fasta file with sequences of all components
        needs Config.py and Trans.py objects
        
    """
    
#@TODO: clean this function; here user should be able to easily convert pdb files into complex or check components correctness!!!

    inp.check_fastanames_duplications()
    
    #if no shape descriptor is provided do not calculate shape fill!!
    if inp.mapfile == None and inp.saxsfile == None:
        con.freespace_penalty = [0.,0.]
        con.density_penalty = [0., 0.]
    index, components_index = 0, {}
    for struct in inp.structures:
        pyrystruct, fasta_seq = inp.check_input(struct)
        #generate subunits objects
        component = trans.get_component(pyrystruct, fasta_seq, con)
        components_index[component.pyrystruct.chain] = index
        index += 1
        
    ##set componenent with no structure as Pyry3D components    
    components_no_struct = inp.get_seq_with_no_struct()
    for seq_no_struct in components_no_struct:
        fasta_seq = seq_no_struct.seq
        seq_name = seq_no_struct.name
        pyrystruct = inp.check_seq_no_struct_input(seq_no_struct)
        trans.get_component_no_struct(seq_no_struct, pyrystruct, con)
        components_index[seq_no_struct.name.split("_")[0]] = index
        index += 1
    #not used yet, will be implemented in next generation of PyRy3D
    con.components_indexes = components_index

    #check whether all chains defined in covalent bonds exist in provided structures
    if len(con.linked_components) != 0: check_corretness_of_covalent_bonds(con)

def initialize():
    """
       checks whether there are any old ouputs in output folder provided by the user;
       if so, it deletes them; otherwise PyRy3D creates new folder to store simulation results
    """
    
    # ----remove decoys from previous PyRy run ---------
    if os.path.exists(str(inp.outname)) == True:
        rmtree(str(inp.outname))
        os.mkdir(str(inp.outname))
    else: os.mkdir(str(inp.outname))
    
    #------ set log, trafl, outfolder names ----------
    logfile.set_filename(inp.outname+"/pyry.log")
    if inp.traflfile != None:
        trafl.set_filename(inp.outname+"/"+inp.traflfile+".trafl")
        
    if inp.movehistory_file:
        movehist.set_filename(inp.outname+"/"+inp.movehistory_file)
        
    outfolder.set_dirname(inp.outname)

def run_simulation(first_complex, simul_params, interactions, density_map, traflfile):
    """
        calls Monte Carlo simulation for complex building
    Parameters:
    -----------
        first_complex       : PyRycomplex object
        interactions        : object of interaction class storing experimental data such as distances, symmetry etc.
        simul_params        : simulation parameters (number of steps,
                              number of out complexes etc.)
        density_map         : map object representing cryo-EM map
        traflfile           : name of traflfile to store trajectory
    Returns:
    ----------
        complexes           : list of PyRyComplexes after simulation (ready models)
    """
#-----------perform simulation!-------------------------------------------------
    #decide which mutations are available for the system
   
    if len(first_complex.movable) == 0: raise InputError("There are no components to mutate. I have northing to do!")
    
    set_available_mutations(first_complex, simul_params)
    
    
    sim = Simul (simul_params)
    sim.setSimulationMethod (simul_params.simmethod)   #@todo Load configuration value
    if simul_params.simmethod == "genetic": #genetic
        sim.setReductionMethod (simul_params.reductmethod)
        sim.setMaximumPoolSize (simul_params.maxpoolsize) #100)	#czemu hard-coded? - Mateusz Susik       
    elif simul_params.simmethod == "simulatedannealing" or \
         simul_params.simmethod == "replicaexchange":
        sim.setTemperature (simul_params.anntemp)
    if simul_params.simmethod == "replicaexchange":
        sim.setReplicaExchangeSteps(simul_params.replica_exchange_freq)
        sim.setReplicaTemperatures(simul_params.replica_temps)
    
    sim.setParameterScalingBoundries(simul_params.param_scaling_ranges)
    
    sim.setParameterScalingValues([ \
    simul_params.param_scaling_range1, simul_params.param_scaling_range2, \
    simul_params.param_scaling_range3 ])
   
    sim.setMutationFrequencies(simul_params.rotation_freq, simul_params.rotation_cov_freq, simul_params.translation_freq, \
                          simul_params.exchange_freq, simul_params.exchangeandsample_freq, simul_params.simul_dd_freq, simul_params.translation_all_freq,\
                          simul_params.rotation_all_freq, simul_params.rotation_whole_freq)
    
    sim.setReheatingParams(simul_params.reheat, simul_params.reheat_rejections)

    sim.setScalingOption (simul_params.param_scaling)
    sim.setStartingComplex (first_complex)
    sim.setIterationCount (simul_params.simul_steps)
    sim.setIndicators(simul_params)
    sim.setResultCount (simul_params.struct_nr)
    sim.setInteractions (interactions)
    sim.setDensityMap (density_map)
    sim.setStepsToSave(simul_params.niter)
    sim.setTraflfile(traflfile)
    sim.setServer(inp.optimized)
    if inp.scoreplot: sim.setScoresPlot()

    best_complex, last_complex = sim.start()    
    
#---------take simulation output according to user's wishes---------------------
    #complexes = sim.getResult()
    #logfile.write_message("Number of rejected complexes is "+str(sim.rejected))
    return best_complex, sim, last_complex

def save_fullatom_bestmodel(best_complex):
    """
    Used to build full atom representation of model obtained by PyRy3D.
    The method can be usefull mainly when a user applied coarse grain representation for model calculation i.e. to speed up the modeling process
    It takes the original structures provided by the user (from input .pdb files) and applies all moves performed
    during simulation to them in order to receive full atom representation for the best complex (with the lowest score)
    Parameters:
       best_complex  : PyRyComplex object with PyRy3D model with the best score (minimal penalty)
    """
    
    print "Saving best model structure in full atom representation. It might take a while. Please be patient...."
    print "The best complex generated by PyRy3D program is ", best_complex.simulation_score, best_complex.iteration_index
    best_complex.save_pdb("best_simulation_complex.pdb")
    trans2 = Trans()
    components_names = []
    for orig_comp in inp.structures:
        pyrystruct, fasta_seq = inp.check_input(orig_comp)
        #generate subunits objects
        trans2.get_component(pyrystruct, fasta_seq, con, "fullatom")
        
    ##set componenent with no structure as Pyry3D components    
    for seq_no_struct in inp.components_no_struct:
        fasta_seq = seq_no_struct.seq
        seq_name = seq_no_struct.name
        pyrystruct = inp.check_seq_no_struct_input(seq_no_struct)
        trans2.get_component_no_struct(seq_no_struct, pyrystruct, con)
        
    original_complex = PyRyComplex(trans2.components)
    original_complex.save_pdb('original_structure_fullatom.pdb')
    
    
    index = 0
    for component in original_complex.components:
        #print "MoveSet history for best scored component:!", best_complex.components[index].pyrystruct.chain, \
        #"rotations", best_complex.components[index].rot_history, "translations" ,best_complex.components[index].trans_history
        
        component.simulate_disorder(component.pyrystruct.struct)
        original_complex.get_alfa_atoms(component)
        original_complex.add_complex_atoms(len(list(component.pyrystruct.struct.get_atoms())))
        original_complex.add_complex_alfa_atoms(len(component.alfa_atoms))
        
        rothist = deepcopy(best_complex.components[index].rot_history)
        translations = deepcopy(best_complex.components[index].trans_history)
        for rot in rothist:
    #@TODO tu dodac rotacje wzgledem wiazania jesli axis==L!!!!
            #if rot[0] == "L":
            #    component.rotate_around_covalent_bond(rot[1], points = rot[2])
            #else:
            if   rot[0][0] != 0.: component.rotate("", rot[0][0], "X")
            elif rot[0][1] != 0.: component.rotate("", rot[0][1], "Y")
            elif rot[0][2] != 0.: component.rotate("", rot[0][2], "Z")
            elif rot[1][0] != 0.: component.rotate_around_covalent_bond(rot[1][0][0], False, rot[1][1][0], rot[1][2][0])
                #component.rotate(rot[1], rot[0])
        for trans in translations:
            print "translate", trans
            if trans[0] != 0 and trans[1] != 0 and trans[2] != 0 :
                component.translate(trans)
        #component.translate(best_complex.components[index].trans_history)
        index += 1
    original_complex.save_pdb(1, inp.fullatom_file)
    
def set_available_mutations(first_complex, simul_params):
    """
       method to check whether all mutations are properly assigned to both
       complex and individual complex components
    Parameters:
       first_complex : PyRyComplex object
       simul_params  : contains parameters of simulation provided by the user
    """
    
    #no exchange when complex consists of one component with unlimited moves
    if len(first_complex.free) < 2:
        #available_mutations.remove("Exchange")
        for index in first_complex.movable:
            if "exchange" in first_complex.components[index].moves.allowed_transform:
                first_complex.components[index].moves.allowed_transform.remove("exchange")
            if "exchangeandsample" in first_complex.components[index].moves.allowed_transform:
                first_complex.components[index].moves.allowed_transform.remove("exchangeandsample")
    for component in first_complex.components:
        if len(simul_params.disabled_mutations) > 0:
            for dis_mut in simul_params.disabled_mutations:
                if dis_mut in component.moves.allowed_transform:
                    component.moves.allowed_transform.remove(dis_mut)
        if len(first_complex.with_disorders) == 0:
            if "SimulateDisorder" in component.moves.allowed_transform:
                component.moves.allowed_transform.remove("SimulateDisorder")
                
def save_movehistory(best_complex, movehist):
        
        for comp in best_complex.components:
            movehist.write_message("COMPONENT\t"+comp.pyrystruct.chain)
            movehist.write_message("TRANSLATIONS\t"+str(comp.trans_history))
            movehist.write_message("ROTATIONS\t"+str(comp.rot_history))
        movehist.write_file()
                                

if __name__=='__main__':
    
    doc = """
    PyRy3D 
    program for macromolecular docking into cryoEM maps
    
    (c) 2010 by Joanna M. Kasprzak
    
    usage: python pyry3d.py --help
    """
    #seed(0)

    print doc

# ---- call main modules objects! -----------------------------------------    
    inp   =    Input()  #represents data delivered by the user
    trans =    Trans()  #represents user data tranformed into Pyry objects
    con   =    Config() #represents user defined simulation parameters

    options = inp.get_data()
    
#-----------set outfolder, logfile and trafl file -----------------------------
    initialize()
    
    warnings.filterwarnings(action='ignore', module='Bio.PDB')
    
#-----------take simulation parameters or decide to use default values---------

    #no file was provided describing complex shape (density map nor saxs ab initio model)
    if inp.mapfile == None and inp.saxsfile == None: shape_descriptor = False
    else:  
        if inp.mapfile: shape_descriptor = "map"
        else: shape_descriptor = "saxs"
        
    con.parse_config_file(inp.config, inp.curvefile, shape_descriptor)
          
    con.set_movehistory_file(inp.movehistory_file)
    
    if inp.movehistory_file:
        con.movehistory = inp.movehistory_file
    
#---geting sequences, structures, restraints and density map provided by the user
    get_input_data(inp, con, trans)
    if options.seq_filename == None and con.identify_disorders == True:
        raise InputError("To identify disordered regions in components you have to provide fasta file with full-length sequences")
   
# ---------- create first complex --------------------------
    first_complex = create_first_complex(trans.components, con, inp.traflfile)
#------------generate interactions between complex components
    density_map = trans.get_map(inp.saxsfile, inp.mapfile, con, first_complex)

    trans.get_interactions(inp.restraints, inp.symrestraints, trans.components, density_map)
    #first_complex.interactions = trans.interactions

#------------locate components randomly inside the density map----------
    calculate_simulbox_statistics(first_complex, density_map, inp, trans, con)
    arrange_components_inside_simulbox(first_complex, density_map, trans, con)
    
#------------perform simulation!-----------------------------------------------
    #import cProfile
    
    if con.simul_steps != 0:
        best_complex, sim, last_complex = run_simulation(first_complex, con, trans.interactions, density_map, inp.traflfile)
       # best_complex, sim, last_complex = cProfile.run('run_simulation(first_complex, con, trans.interactions, density_map, inp.traflfile)')
    
    if inp.fullatom_file: save_fullatom_bestmodel(best_complex)
    
# ---------- save log messages -------------------------------------------     
    logfile.write_file()
    #if inp.movehistory_file: save_movehistory(best_complex, movehist) #movehist.write_file()
    if inp.movehistory_file: save_movehistory(last_complex, movehist) #movehist.write_file()
    
    if inp.scoreplot: draw_scores_plot(inp.scoreplot, sim)
    
    
    
