#!/usr/bin/env python
from __future__ import division
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

from copy                                  import deepcopy
#numpy
from numpy                                 import array, dot, arange, append
#math
from math                                  import sin, cos, radians, sqrt, fabs, degrees
#random
from random                                import randint, sample
#biopython
from Bio                                   import PDB, pairwise2
from Bio.PDB                               import Entity
from Bio.PDB.Atom                          import Atom 
from Bio.PDB.Residue                       import Residue
from Bio.PDB.Chain                         import Chain
from Bio.PDB.Model                         import Model
from Bio.PDB.Structure                     import Structure
#PyRy
from Modules.Constans.Constans             import ATOM_RADII, MOLWEIGHTS, \
                                            AA_WEIGHTS, DNA_WEIGHTS, RNA_WEIGHTS, LIGANDS
from Modules.Error.Errors                  import InputError, ComponentError, Component_movesError
from Modules.Input.PyRyStructure           import PyRyStructure, PyryAtom
from Modules.Trans.ComponentRepresentation import ComponentRepresentation
from Modules.Trans.VolumeSimulator         import Grapes, Disordered_Fragment, point
from Modules.Constans.Logfile              import logfile


"""@TODO:
check names of methods and variables
remove code redundancy (has_changed, check_if_changed)
check private vs public methods!!
__str__ method
missing docstrings, descriptions etc
"""


class Component(object):
    """
        class that stores information about all complex components,
        which are:
            chains with assigned structure
            sequences with calculated volume
    """
              
    def __init__(self):
        self.pyrystruct      = None          #PyRyStructure object
        self.fasta_seq       = None          #Bio object of fasta sequence
        self.pyryatoms       = []            #list of pyryatoms objects
        self.volume          = 0.0           #vomume of given component
        self.vollist         = []            #to store MW of residues
        self.mass_centre     = [0.,0.,0.]    #stores x,y,z coordinates of centre of mass
        self.rot_history     = [[[0.,0.,0.],[0]]]            #all rotations during simulation
        self.rot_history_sum = {"X":0.0, "Y": 0.0, "Z": 0.0, "L":0.0}        #sum of all rotations for given component
        self.trans_history   = [[0.0, 0.0, 0.0]]  #all translations -||-
        self.moves           = None
        self.alfa_atoms      = []             #stores calfa atoms
        self.covalent_bonds  = []             #stores indexes of components bound with given component by covalent bond
        
        #--------- simulation parameters for given component only --------------        
        self.outbox_atoms   = 0               #number of atoms outside a simulbox
        self.mapcells_nr     = 0.              #number of mapcells neccsery to describe particular component (default, constant value)
        self.taken_mapcells = {}
        
        
        self.taken_mapcells_coords = []
        self.taken_density  = 0.              #sum of density values occupied by a particular component
        self.outmap_cells   = 0.              # -||- but in given simulation step
        self.collided_atoms = []              #list of atoms collided in given simul iteration
        self.outbox_atoms_li = []             #list of atoms outside simul box
        self.inmap_atoms     = []             #list of atoms numbers inside the map
        
        self.sa_interactions_score = 0.       #temp score for surface access restraint, later summed up with other restraints
        self.pd_interactions_score = 0.       #temp score for poind distance restraint, later summed up with other restraints
       
        #---------information about disorder -----------------------------------
        self.disorders      = []                 #list of disordered regions
        self.is_positioned  = False            #just to check whether a component has been positioned randomly in the map or not
        #---------------------------------------------------------------------------

        self.is_outmap      = False          #True if a Component (all its atoms) are outside the density map
        self.is_outbox      = False          #True if a Component (all its atoms) are outside the simulation box

        self.repr_type      = None          #ca, cacb, 3p, fa, sphere, ellipsoid
        self.sphere_mw      = 0.            #temporal value stored only for components with sphere representation
	self.sphere_radius  = 0.            #temporal value stored only for components with sphere representation
	
    def __str__(self):
        return "%s " % (self.pyrystruct)
        
    def add_alfa_atom(self, atom):
        self.alfa_atoms.append(atom) #pyryalfa)
        
    def add_component_moveset(self, simul_params):
        """
        """
        self.moves = Component_moves()
        movable = simul_params.movable
        found = False
        for comp in movable:
            if comp.comp_name == self.pyrystruct.chain:
                if comp.state == "fixed":
                    found = True
                    self.moves.set_fixed_component()
                    logfile.write_message("MOVEMENTS of COMPONENT "+comp.comp_name+" "+comp.state)
                    break
                elif comp.state == "movable":
                    found = True
                    self.moves.set_movable_component(comp, simul_params)
                    logfile.write_message("MOVEMENTS of COMPONENT "+comp.comp_name+" "\
                                         +comp.state+str(comp.max_trans_vector)+str(comp.xmax_rot)+\
                                         " "+str(comp.ymax_rot)+" "+str(comp.zmax_rot))
                    break   

        if found == False:
            self.moves.set_default_moves(simul_params)
        
    def add_covalent_bonds(self, con, pyrycomplex):
        """
        adds information about connections with other complex components by covalent bonds
        connected components are mutated together
        """
        
        #print "--", self.pyrystruct.chain, con.components_indexes, con.linked_components
        #independent component with no covalent bonds, can move independently
        if self.pyrystruct.chain not in con.linked_components.keys():
                self.covalent_bonds = []
        else:
            for covbond in con.linked_components[self.pyrystruct.chain]:
                #check correctness of covalent bonds atoms defined by the user
                if self.pyrystruct.struct[0][self.pyrystruct.chain].has_id((' ', covbond.atom1[0], ' ')):
		    res1 = self.pyrystruct.struct[0][self.pyrystruct.chain][covbond.atom1[0]]
                else:
                    raise InputError("Component %s does not have residue %s"%(self.pyrystruct.chain, covbond.atom1))
                #
                if self.pyrystruct.struct[0][self.pyrystruct.chain].has_id((' ', covbond.atom2[0], ' ')):
                    res2 = self.pyrystruct.struct[0][self.pyrystruct.chain][covbond.atom2[0]]
                else:
                    raise InputError("Component %s does not have residue %s"%(self.pyrystruct.chain, covbond.atom2))
                
                if res2.has_id(covbond.atom1[1]): pass
                else: raise InputError("Residue %s does not have atom %s"%(covbond.atom1[0], covbond.atom1[1]))
                                  
                if res2.has_id(covbond.atom1[1]): pass
                else: raise InputError("Residue %s does not have atom %s"%(covbond.atom2[0], covbond.atom2[1]))
                
                self.covalent_bonds.append(covbond)
                for chain in covbond.chains:
                    covbond.chains_indexes.append(pyrycomplex.get_component_index_by_chain(chain))                #(con.components_indexes[chain])
            #print "@@@@@", self.covalent_bonds, con.components_indexes, covbond.chains_indexes #, type(covalent_bonds)
         
    def add_disordered_region(self, disorder):
        """
        adds Disordered_Fragment instance to list of disordered fragments for
        a particular component
        """
        self.disorders.append(disorder)
        
    def build_rotation_matrix(self,point1,point2, angle):
        """
        defines a rotation matrix for toration around line defined by
        points point1 and point2 by angle
        
        source:
        http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/
        """
        u = (point2[0]-point1[0])
        v = (point2[1]-point1[1])
        w = (point2[2]-point1[2])
        
        if u == 0. and v ==0. and w == 0:
            raise InputError("To define line you need two different points in 3D space. Please change coordinated for \
                             COVALENT_BOND definition in configuration file")
        
        uu = u*u
        vv = v*v
        ww = w*w
                
        a,b,c = point1[0], point1[1], point1[2]
        
        cosT = cos(angle)
        sinT = sin(angle)
        oneMinusCosT = 1-cosT
        
        l2 = uu + vv + ww
        l = sqrt(l2)
        
        #Build the matrix entries element by element.
        m1a = (uu + ((vv + ww) * cosT))/l2
        m1b = ((u*v * oneMinusCosT - w*l*sinT))/l2
        m1c = ((u*w * oneMinusCosT + v*l*sinT))/l2
        m1d = (((a*(vv + ww) - u*(b*v + c*w)) * oneMinusCosT\
            + (b*w - c*v)*l*sinT))/l2
        
        m1 = [m1a, m1b, m1c, m1d]
        
        m2a = (u*v * oneMinusCosT + w*l*sinT)/l2
        m2b = (vv + ((uu + ww) * cosT))/l2
        m2c = ((v*w * oneMinusCosT - u*l*sinT))/l2
        m2d = (((b*(uu + ww) - v*(a*u + c*w)) * oneMinusCosT\
            + (c*u - a*w)*l*sinT))/l2
        m2 = [m2a, m2b, m2c, m2d]
        
        m3a = (u*w * oneMinusCosT - v*l*sinT)/l2
        m3b = (v*w * oneMinusCosT + u*l*sinT)/l2
        m3c = (ww + ((uu + vv) * cosT))/l2
        m3d = (((c*(uu + vv) - w*(a*u + b*v)) * oneMinusCosT\
            + (a*v - b*u)*l*sinT))/l2
        m3 = [m3a, m3b, m3c, m3d]
                
        return m1,m2,m3
        
    def calculate_component_volume(self):
        """
        assigns volume value to given component. the volume is based on sequence
        molecular masses only!
        """
#----include in tests!! -------------------      
        if self.pyrystruct.moltype.lower()   == 'protein':
            self.__calc_protein_volume()
        elif self.pyrystruct.moltype.upper() == 'DNA':
            self.__calc_nucl_volume('dna')
        elif self.pyrystruct.moltype.upper() == 'RNA': 
            self.__calc_nucl_volume('rna')
            
    def check_rotation_limits(self, axis):
        """
        method designed for components with movements limited by the user; it checks
        what rotations were applied for the components and changes ranges for furher
        movements in order to prevent from exceeding user's limitations
        """
        #for rotations around line limitations are as for X axis
        #print "Checking rotation limits", self.moves.rot_ranges[axis][0], self.moves.rot_ranges_def[axis][0], self.rot_history_sum[axis]
        if self.moves.rot_ranges_def[axis][0]  >=  self.rot_history_sum[axis]:
            self.moves.rot_ranges_sum[axis][0] = 0
        elif self.moves.rot_ranges_def[axis][1]  <=  self.rot_history_sum[axis]:
            self.moves.rot_ranges_sum[axis][1] = 0
        else:
            i1, i2 = self.__return_axis_index(axis)
            self.moves.rot_ranges_sum[axis][1] = self.moves.rot_ranges_sum[axis][1] - self.rot_history[-1][i1][i2]
            self.moves.rot_ranges_sum[axis][0] = -(self.rot_history[-1][i1][i2]+fabs(self.moves.rot_ranges_sum[axis][0]))  
            if round(fabs(self.moves.rot_ranges_sum[axis][1]),5) + round(fabs(self.moves.rot_ranges_sum[axis][0]),5) > round(fabs(2* self.moves.rot_ranges_def[axis][0]),5):
                raise Component_movesError("This should never happen while rotating!")
        #print "--Checking rotation limits--", self.moves.rot_ranges_sum, self.rot_history_sum[axis], axis, self.pyrystruct.chain
        
    def __return_axis_index(self, axis):
        if   axis == "X": return 0, 0
        elif axis == "Y": return 0, 1
        elif axis == "Z": return 0, 2
        elif axis == "L": return 1, 0
        
    def check_translation_limits(self):
        """
        method designed for components with movements limited by the user; it checks
        what translations were applied for the components and changes ranges for furher
        movements in order to prevent from exceeding user's limitations
        """
        #print "Checking transl limits", self.moves.trans_ranges, self.trans_history, self.moves.trans_ranges_def
        for index in xrange(0,3):
            trans_hist_index = self.summarize_translations(index)
            if self.moves.trans_ranges_def[index][0] >=  trans_hist_index:
                self.moves.trans_ranges_sum[index][0] =  0 
            elif self.moves.trans_ranges_def[index][1] <=  trans_hist_index: 
                self.moves.trans_ranges_sum[index][1] = 0
            else:
                self.moves.trans_ranges_sum[index][1] = self.moves.trans_ranges_sum[index][1] - self.trans_history[-1][index] #self.trans_history[index]
                self.moves.trans_ranges_sum[index][0] = -(self.trans_history[-1][index] + fabs(self.moves.trans_ranges_sum[index][0]))
                if round(fabs(self.moves.trans_ranges_sum[index][1]),5) + round(fabs(self.moves.trans_ranges_sum[index][0]),5) > round(fabs(2* self.moves.trans_ranges_def[index][0]),5):
                    raise Component_movesError("This should never happen while translating!")
        #print "---Checking transl limits---", self.pyrystruct.chain, self.moves.trans_ranges, self.moves.trans_ranges_sum, self.trans_history[-1]
        
    def choose_rotation_points(self):
        """
        assigns atoms declared as covalent bonds as points defining rotationvector
        
        if a given component has more than one bond, a method selects one bond randomly
        """
        if len(self.covalent_bonds) == 0:
            at1 = list(self.pyrystruct.struct.get_atoms())[0]
            at2 = list(self.pyrystruct.struct.get_atoms())[-1]
            return at1, at2, None
        else:
            random_bond_index = randint(0, len(self.covalent_bonds)-1)
            random_bond = self.covalent_bonds[random_bond_index]
            
            res1 = list(self.pyrystruct.struct.get_residues())[random_bond.atom1[0]]
            res2 = list(self.pyrystruct.struct.get_residues())[random_bond.atom2[0]]
            
            at1 = res1[random_bond.atom1[1]]
            at2 = res2[random_bond.atom2[1]]                  

            #print "selected atom", at1.serial_number, at2.serial_number, res1.id, res2.id
            return at1, at2, random_bond_index
        
    def summarize_translations(self, index):
        """
        sumarizes all translations in given direction defined by index:
        0 - "X", 1- "Y", 2 - "Z"
        """
        trans_sum = 0.
        for trans in self.trans_history:
            trans_sum += trans[index]
        return trans_sum
        
    def __choose_simulated_region(self, option):
        """
        determines which regions will be mutated in a particular simulation step
        during simulation one region is randomly selected; during complex
        formation all regions must be simulated
        """
        #random selection of mutated region
        if option == "simul": 
            mutate_regions = []
            counter = randint(1, 2**len(self.disorders) - 1)
            index = 0
            while counter:
                if counter % 2:
                    mutate_regions.append(index)
                ++index
                counter //= 2
        else:
            mutate_regions = xrange(0, len(self.disorders))
        return mutate_regions
        
    def clean_collided_atoms(self):
        self.collided_atoms = []
        
    def clean_taken_mapcells(self):
        """
        removes all elements from taken_mapcells attribute method called before
        new coordinates in map are calculated
        """
        self.taken_mapcells = {}
        self.taken_mapcells_coords = []
        
    def clear_pyryatoms(self):
        """
           remove a pyryatoms list to free space and computational time
           during simulation
        """
        del self.pyryatoms

    def find_disordered_fragments(self):
        """
        check if structure sequence is a part of fasta sequence and identifies
        sequence fragments with no structure
        
        uses Bio.pairwise2
        """
        #print "&&&", self.pyrystruct.sequence
        try: aligned = pairwise2.align.globalxs(self.fasta_seq.seq, self.pyrystruct.sequence, \
                                           -0.1, 0, penalize_end_gaps=False)
        except:
            print "Your component is to large to identify disorders!!!"
            return []
        gaps = self.refind("-", aligned[0][1])
        #print "FASTA\t",self.fasta_seq.seq
        #print "STRUCT\t", self.pyrystruct.sequence
        dd_frag_positions = self.__identify_dd_positions(gaps)
        
        if len(dd_frag_positions) != 0:
            self.__find_if_dd_in_struct(dd_frag_positions)
        dd_frag_seqs = self.__get_dd_frag_seqs(dd_frag_positions)
        
        return aligned
    
    def __find_if_dd_in_struct(self, dd_frag_positions):
        
        """
        [[1,2,3], [5,6,7], [10,11]]
        """
        resi_nrs = []
        for resi in self.pyrystruct.struct.get_residues():
            resi_nr = resi.id[1]
            resi_nrs.append(resi_nr)
	    resi_nrs.sort()
	#print resi_nrs
                
        dd_indexes = []
        for frag in dd_frag_positions:
            for f in frag:
                dd_indexes.append(f)
	#print "^^^", dd_indexes
        
        #print resi_nrs, len(list(self.pyrystruct.struct.get_residues())), len(self.pyrystruct.sequence)
        #print "DDfrag", dd_frag_positions
        #print "Checking", dd_indexes, resi_nrs, len(list(self.pyrystruct.struct.get_residues())), len(resi_nrs)
        for i in xrange(1, len(self.fasta_seq)):
            #if disordered fragment is too long
            if (i in resi_nrs) and (i in dd_indexes):
                raise InputError(str(i)+"residue is duplicated in fasta and pdb sequence. \n FASTA seq:"+self.fasta_seq.seq+"\n SEQ in PDB file"+self.pyrystruct.sequence)
            #if disordered fragment is too short
            if (i not in resi_nrs) and (i not in dd_indexes):
                raise InputError(str(i)+"residue is missing when comparing fasta and pdb sequence in chain "+(self.pyrystruct.chain)+"\n FASTA seq:"+self.fasta_seq.seq+"\n SEQ in PDB file"+self.pyrystruct.sequence)
    
    def refind(self, substring, seq, start=0):
        """
        Returns all start positions of substring in seq.
        by Kristian Rother
        Returns:
        --------
        a list with indexes of all gaps in fasta sequence
        e.g. for:
        ---joanna--maria-kasprzak-
        will return
        [0,1,2,9,10,16,25]
        """     
        res = []
        for i in range (0, len(seq)):
            if (seq[i] == substring):  res.append(i)
        return res
            
    def __identify_dd_positions(self, gaps):
        """
        identifies all unaligned regions and assignes them to separate lists of
        indexes
        Parameters:
        ----------
        gaps - lists of indexes of all gaps in aligment
        Returns:
        --------
        a list with indexes of all unaligned fragments in fasta sequence
        e.g. for:
        ---joanna--maria-kasprzak-
        will return
        [[0,1,2], [9,10], [16], [25]]
        """
        if len(gaps) != 0:
            sequence, results = [], []
            sequence.append(gaps[0]+1)
            for index in range(0, len(gaps)-1):
                if gaps[index] +1 == gaps[index+1]:
                    sequence.append(gaps[index+1]+1)
                else:
                    results.append(sequence)
                    sequence = []
                    sequence.append(gaps[index+1]+1)
            results.append(sequence)
            return results
        else: return []
    
    def __get_dd_frag_seqs(self, dd_frag_positions):
        """
        Having a list with all disordered fragments location a method retrieves
        a start and end position of each fragment,
        retrieves corresponding fasta sequence fragment and calls
        Disordered_Fragment object
        """
        for ddfrag in dd_frag_positions:
            start_pos = ddfrag[0]
            stop_pos = ddfrag[-1]
            frag_seq = self.fasta_seq.seq[start_pos-1: stop_pos]

            dd_frag = Disordered_Fragment(start_pos, stop_pos, frag_seq)
            dd_frag.create_simulated_fragment(self.pyrystruct, self.fasta_seq.seq)
            self.add_disordered_region(dd_frag)
            print "Disorder found!!", dd_frag
            
    def rotate(self, param, rangle, axis="X", history = True, vector = [0.,0.,0.]):
        """
        performs rotation and changes structure coordinates according to rotation angle
        The rotaion matrix this - z axis rotation when looking towards the origin.  
        Parameters:
        -----------
            angle       :   value from -360 to 360
            axis        :   x, y or z
            history     :   True indicates that method should save this move
        Raises:
        -------
            Cmplx_ComponentsError if wrong axis name or angle value is given
        """

        angle = radians(rangle) #Converse gamma to radians
        if history == True:            
            if param != "rotate_whole":
                    if   axis == "X": self.rot_history.append([[rangle, 0,0],[0]])
                    elif axis == "Y": self.rot_history.append([[0, rangle, 0],[0]])
                    elif axis == "Z": self.rot_history.append([[0,0,rangle],[0]])

            self.rot_history_sum[axis] += angle
            #print "rothist", self.rot_history, len(self.rot_history)
        
        minus_vector = vector if param == "rotate_whole" else -self.mass_centre
        #if param == "rotate_whole":
        #        print "Rotate whole", self.pyrystruct.chain, degrees(angle), axis
        #else:
        #print "Rotate", self.pyrystruct.chain, degrees(angle), axis, minus_vector #, self.rot_history


        
        if   axis == 'X':  rot = array([[1,0,0],[0, cos(angle), sin(angle)],\
                                                [0,-sin(angle),cos(angle)]])
        elif axis == 'Y':  rot = array([[cos(angle), 0, -sin(angle)],[0,1,0], \
                                        [sin(angle), 0, cos(angle)]])  
        elif axis == 'Z':  rot = array([[cos(angle), sin(angle), 0],\
                                        [-sin(angle), cos(angle), 0], [0,0,1]])
        else: raise ComponentError("No such axis %s !!"%(axis))
        
        for atom in self.pyrystruct.struct.get_atoms():
            rotated_coordinates = dot((atom.coord+minus_vector),rot)+(-minus_vector)
            atom.coord = rotated_coordinates
            
    def rotate_around_covalent_bond(self, rangle, history = True, points=[], cov_index = None):
        """
        rotates a component around a line defined by two points.
        """
        angle = radians(rangle) #Converse gamma to radians
        if points:
            point1 = points[0]
            point2 = points[1]
        else:
            point1, point2, cov_index = self.choose_rotation_points()
        
        if history == True:
            self.rot_history.append([[0,0,0],[[rangle], [point1.get_serial_number()], [point2.get_serial_number()]]])
            self.rot_history_sum["L"] += angle
                
        m1,m2,m3 = self.build_rotation_matrix(point1.coord,point2.coord, angle)
        
        for atom in self.pyrystruct.struct.get_atoms():
            rotated_coordinates = self.calc_coordinates2(atom.coord,m1,m2,m3)
            atom.coord = rotated_coordinates        
        self.set_center()    
        #print "RotateCOVBOND", self.pyrystruct.chain, degrees(angle), point1.get_serial_number(), point2.get_serial_number(),cov_index #, self.rot_history
        
        return [point1, point2], cov_index
    
    def calc_coordinates2(self, coords,m1, m2,m3):
        """
        """
        x,y,z = coords[0], coords[1], coords[2]
        m11, m12, m13, m14 = m1[0],m1[1],m1[2],m1[3]
        m21, m22, m23, m24 = m2[0],m2[1],m2[2],m2[3]
        m31, m32, m33, m34 = m3[0],m3[1],m3[2],m3[3]
        
        p = m11*x + m12*y + m13*z + m14
        q = m21*x + m22*y + m23*z + m24
        r = m31*x + m32*y + m33*z + m34
    
        new_coords = array([p,q,r])
        return new_coords
    
    def set_center(self):
        """
        for given complex component method calculates estimated center
        which is essential for all geometrical operations like rotation,
        translation etc  
        Returns:
        --------
            [x,y,z] coordinates of the centre
        """
        self.mass_centre = self.pyrystruct.calculate_centre_of_mass()
    
    def set_component(self, pyrystructure, fasta_seq, simul_params, option = None, radius = None):
        """
        assigns PyRyAtoms to structure
        Parameters:
        -----------
            pyrystrusture       : PyRyStructure object - single complex component
            grid_type           : 
                                  'coarse' - coarse grain model
                                  'fullatom' - full atom model
                                  'grid_cubic' - cubic grid
                                  'grid_diament' - diament grid
            radius              : None for non-grid
                                  float for grids, if no value is given, radius is 3A         
        """
        
#TODO: if crude representation: set it as pyrystruct!!
        self.pyrystruct = pyrystructure
                
        self.fasta_seq = fasta_seq
        
        #TODDO go before self.pyrystruct
        self.add_component_moveset(simul_params)
        
        if option: pass
        elif option == None:
            self.__create_component_representation(simul_params.representation) #, radius)
        
        #get disordered fragments!
        if simul_params.identify_disorders == True:
            self.find_disordered_fragments()
        
        
        self.__collect_pyryatoms(simul_params.representation)
        
        atoms = list(self.pyrystruct.struct.get_atoms())
        for at in atoms:
            atom_parent = at.get_parent()
            atom_parent.detach_child(at.get_name())
            for pa in self.pyryatoms:
                if pa.serial_number == at.serial_number:
                    atom_parent.add(pa)
            
        self.set_center()
        self.calculate_component_volume()
        
    def set_component_no_struct(self, fasta_seq, pyrystruct, simul_params):
        """
        assigns PyRyAtoms to structure
        Parameters:
        -----------
            pyrystrusture       : PyRyStructure object - single complex component
            simul_params   
        """
        self.fasta_seq = fasta_seq
        
        volume = Disordered_Fragment(1, len(self.fasta_seq), self.fasta_seq)
        volume.create_simulated_volume(1, self.fasta_seq, pyrystruct)
        
        g = Grapes()
        #self.mass_centre = srodek mapy!!!!
        g.set_volume_simulation_parameters(volume, self, pyrystruct.struct, pyrystruct.moltype, self.mass_centre)
        
        res = g.generate()        
        volume.add_pseudoatoms_to_structure(res, pyrystruct.moltype)
        
        
        self.pyrystruct = pyrystruct
        self.pyrystruct.set_structure(volume.fragment_lattice)  
        self.add_component_moveset(simul_params)
        
        #---------get disordered fragments!----------------
        self.add_disordered_region(volume)
        print "disordered component with no structure found!!", volume #disorder

        self.set_center()
        self.calculate_component_volume()
	
        
    def set_mapcells_nr(self, filled_cells):
        """
        sets how many simgrid cells are taken by particular complex component
        in given simulation step
        """
        #self.filled_cells = filled_cells
        self.mapcells_nr = filled_cells
        
    def set_gridcells(self, gridcells):
        """
        sets number of grid cells occupied by a component
        """
        self.gridcells_nr = gridcells

    def simulate_disorder(self, struct = False, option = None):
        """
        performs simulate disorder mutation. first randomly chooses fragment to be
        changed than runs VolumeSimulator in chosen mode:
        -- grape like (for regions at terminis of proteins and for sequences with no structure)
        -- chain like (for disordered fragments inside a molecule and for nucleic acids)
        """
        
        #print "Sim DD for component", self.pyrystruct.chain
        if len(self.disorders) == 0: return 0    #no disordred regions in this component   
        mutate_regions = self.__choose_simulated_region(option)

    #return ddstruct
        self.__simulate_volume(struct, mutate_regions)
        return 1
    
    def set_pd_restraint_score(self, pd_rest_score):
        self.pd_interactions_score = pd_rest_score
    
    def set_sa_restraint_score(self, sa_rest_score):
        self.sa_interactions_score = sa_rest_score
    
    def __simulate_volume(self, struct, mutate_regions):
        """
        """
        modified_dd_frags = []
        
        for index in mutate_regions:
            dd_frag = self.disorders[index]
            dd_frag.set_modeling_disordered_fragment(self, struct, self.pyrystruct.chain) 

            g = Grapes()
            g.set_volume_simulation_parameters(dd_frag, self, struct, self.pyrystruct.moltype, self.mass_centre) 

            res = g.generate(struct) #self.pyrystruct.struct)
            dd_frag.add_pseudoatoms_to_structure(res, self.pyrystruct.moltype)
            
#@TODO might be VERY error prone!!
            #print "simulate volume", dd_frag.fragment_type, dd_frag.start_pos, dd_frag.stop_pos+1, self.pyrystruct.chain
            if dd_frag.fragment_type != "simulated_volume":
                #print "-----------------", len(res), len(dd_frag.fragment_lattice) #, len(dd_frag.structure)
                dd_frag.add_fragment_to_original_structure(self, struct, dd_frag.start_pos, dd_frag.fragment_type)   #self.pyrystruct.(start_pos)
                #----volume simulation always change center of mass-----
                self.set_center()
            else:
                self.pyrystruct.struct = dd_frag.fragment_lattice
                struct = dd_frag.fragment_lattice
                #----volume simulation always change center of mass-----
                self.set_center()
                
            dd_frag.clean_fragment()
            g.clean_grapes()
        return struct
            
        #print "!!!", dd_frag.sequence, dd_frag.fragment_lattice, len(dd_frag.pseudoresidues)
        
        
    def translate(self, trans_vector, history = True):
        """
        performs translation of given structure via defined vector 
        Parameters:
        -----------
            trans_vector    :   translation vector e.g. [x,y,z]
            history         :   True indicates that method should save this move    
        Returns:
        ---------
            structure coordinates after translation
        """
        self.mass_centre += array(trans_vector)
        
        if history == True:
            self.trans_history.append(list(trans_vector))
            #print "history----", self.pyrystruct.chain, self.trans_history, len(self.trans_history)
        #print "Translate", self.pyrystruct.chain, trans_vector #, self.trans_history
        for at in self.pyrystruct.struct.get_atoms():
            at.coord = array(at.coord) + array(trans_vector)
        self.set_center()
        
    #------Private methods----------
    def __calc_protein_volume(self):
        """
        calculates volume for proteins
        
        (0.73 cm3/g x 10**24A/cm3 x molecular weight g/mole)/(6.02 x 10**23 molecules/mole)
        
        and results in a protein volume of approximately:
        (1.21 x MW) A3/molecule
        """
        #for resi in self.pyrystruct.sequence:
        for resi in self.fasta_seq.seq:
            resi_weight =  AA_WEIGHTS[resi.upper()]
            self.volume += resi_weight
            self.vollist.append(resi_weight)
        self.volume =  self.volume
        
    def __calc_nucl_volume(self, type):
        """
        Calculates volume for RNA and DNA
        Parameters:
        -----------
            type: 'rna' or 'dna'
        """
        if type == "dna":
            #for resi in self.pyrystruct.sequence:
            for resi in self.fasta_seq.seq:
                resi_weight = DNA_WEIGHTS[resi.upper()]
                self.volume += resi_weight
                self.vollist.append(resi_weight)
            self.volume = self.volume - 61.96
            self.vollist.append(-61.96)
        elif type == "rna":
            #for resi in self.pyrystruct.sequence:
            for resi in self.fasta_seq.seq:
                resi_weight = RNA_WEIGHTS[resi.upper()]
                self.volume += resi_weight
                self.vollist.append(resi_weight)
            self.volume = self.volume + 159.00
            self.vollist.append(159.00)
    
    def __collect_pyryatoms(self, repr_type):
        """
        collects all PyRyAtoms and sets its attributes
        """
        for at in self.pyrystruct.struct.get_atoms():
            hetname = at.get_parent().id[0].split("_")
            if hetname[0] == "H":
                if hetname[1].strip().upper() in LIGANDS: continue 
            new_atom = PyryAtom(at.name, at.coord, at.bfactor, at.occupancy, \
                                    at.altloc, at.fullname, at.serial_number)
            if repr_type != "sphere":
		new_atom.assign_vdw()
                new_atom.assign_molweight()
	    else:
		new_atom.vdw = self.sphere_radius
		new_atom.molweight = self.sphere_mw
		del self.sphere_mw
		del self.sphere_radius
		
            self.pyryatoms.append(new_atom)
            self.pyryatoms.sort(key=lambda x:x.serial_number)
            
    def __create_component_representation(self, repr_type): #, radius = None):
        """

        Parameters:
        -----------
        Raises:
        -------
            ComponentError if wrong representation name is provided
        """
        repr = ComponentRepresentation(self.pyrystruct, repr_type)
        if repr_type == 'fa':
	    self.repr_type = "fa"
        elif repr_type == 'ca':
	    self.repr_type = "ca"
            repr = repr.create_coarse_repr(repr_type)
            self.pyrystruct.struct = []
            self.pyrystruct.struct = repr.cg_struct
        elif repr_type == 'cacb':
	    self.repr_type = "cacb"
            repr = repr.create_coarse_repr(repr_type)
            self.pyrystruct.struct = []
            self.pyrystruct.struct = repr.cg_struct
        elif repr_type == '3p':
	    self.repr_type = "3p"
            repr = repr.create_coarse_repr(repr_type)
            self.pyrystruct.struct = []
            self.pyrystruct.struct = repr.cg_struct
	elif repr_type == 'sphere':
	    self.repr_type = "sphere"
            repr  = repr.create_coarse_repr(repr_type)
            self.pyrystruct.struct = []
            self.sphere_radius = repr.radius
            self.sphere_mw     = repr.molmass
            self.pyrystruct.struct = repr.cg_struct
	elif repr_type == 'ellipsoid':
	    self.repr_type = "ellipsoid"
            repr = repr.create_coarse_repr(repr_type)
            self.pyrystruct.struct = []
            self.pyrystruct.struct= repr.cg_struct
        else: raise ComponentError("Sorry, we do not have %s representation available"%(repr_type))
    
class Component_moves(object):
    """
        to store information about allowed moves performed on particular
        complex components from selection:
        rotation/ translation
        rotation_angle
        min_rotation angle, max_rotation_angle
        min_translation_vector, max_translation_vector      
    """
        
    def __init__(self):
        self.state = ''                               #'movable', 'fixed'
        self.allowed_transform = []                   #"rotation", "translation"
        self.rot_axis = []                            #selection of allowed rotation axis eg [x,y], z not allowed
        
        self.rot_ranges = {}                          #{"X": [-10, 10], "Y": [-20,20], "Z": [-30,30]}
        self.trans_ranges = []                        #eg [[-1,10], [-10,1], [-10,5]]
        self.rot_ranges_sum = {}
        self.trans_ranges_sum = []
        
        self.limited = False                          #if component has limited moves
               
    def __str__(self):
        return "%s %s %s %s %s" % (self.state, self.allowed_transform, self.rot_axis, \
                                 self.trans_ranges, self.rot_ranges  )
    
    def check_allowed_transformations(self, move_set):
        """
        """
        if move_set.max_trans_vector[0] != 0 or move_set.max_trans_vector[1] != 0 or move_set.max_trans_vector[2] != 0:
            self.set_transformations("translation")
            self.set_transformations("translation_all")
        if move_set.xmax_rot != 0. or move_set.ymax_rot != 0. or move_set.zmax_rot != 0.:
            self.set_transformations("rotation")
            self.set_transformations("rotation_cov")
            self.set_transformations("rotation_all")
            self.set_transformations("rotation_whole")
            
    def set_default_moves(self, config):
        """
        here assumption is that each component can move with no limitations
        different than those for all components maxrot, maxtrans
        """
        self.state = "movable"
        self.allowed_transform = ["rotation", "rotation_cov", "translation", \
                        "exchange", "exchangeandsample", "SimulateDisorder", "translation_all", \
                        "rotation_all", "rotation_whole"]
        self.rot_axis = ["X", "Y", "Z"]
        
        self.__set_default_rot_ranges(config)
        self.set_translation_vec(config.max_trans_vec, ['10000','10000', '10000'],config, True)
        
    def __set_default_rot_ranges(self, config):
        self.rot_ranges["X"] = [-config.max_rot_angle, config.max_rot_angle]
        self.rot_ranges["Y"] = [-config.max_rot_angle, config.max_rot_angle]
        self.rot_ranges["Z"] = [-config.max_rot_angle, config.max_rot_angle]
        self.rot_ranges["L"] = [-config.max_rot_angle, config.max_rot_angle]
        
        self.rot_ranges_sum["X"] = [-360, 360]
        self.rot_ranges_sum["Y"] = [-360, 360]
        self.rot_ranges_sum["Z"] = [-360, 360]
        self.rot_ranges_sum["L"] = [-360, 360]

        
        self.rot_ranges_def = self.rot_ranges
        #print "ROT ranges NEW!!", self.rot_ranges, self.rot_ranges_def
        
    def set_fixed_component(self):
        """
        """
        self.set_state("fixed")
        
    def set_limited_component(self):
        self.limited = True
        
    def set_movable_component(self, move_set, simparams):
        """
        """
        #comp.set_component_moves()
        self.set_limited_component()
        self.check_allowed_transformations(move_set)
        
        self.set_rotations(move_set)
        
        self.set_transformations("SimulateDisorder")
        move_set.xmax_rot, move_set.xmax_rot_sum = self.set_range_values(move_set.xmax_rot, move_set.xmax_rot_sum, "X", simparams)
        move_set.ymax_rot, move_set.ymax_rotsum = self.set_range_values(move_set.ymax_rot, move_set.ymax_rot_sum, "Y", simparams)
        move_set.zmax_rot, move_set.zmax_rot_sum = self.set_range_values(move_set.zmax_rot, move_set.zmax_rot_sum, "Z", simparams)
        move_set.lmax_rot, move_set.lmax_rot_sum = self.set_range_values(move_set.lmax_rot, move_set.lmax_rot_sum, "L", simparams)
        
        self.rot_ranges["X"] = [-float(move_set.xmax_rot), float(move_set.xmax_rot)]
        self.rot_ranges_sum["X"] = [-float(move_set.xmax_rot_sum), float(move_set.xmax_rot_sum)]
        
        self.rot_ranges["Y"] = [-float(move_set.ymax_rot), float(move_set.ymax_rot)]
        self.rot_ranges_sum["Y"] = [-float(move_set.ymax_rot_sum), float(move_set.ymax_rot_sum)]
        
        self.rot_ranges["Z"] = [-float(move_set.zmax_rot), float(move_set.zmax_rot)]
        self.rot_ranges_sum["Z"] = [-float(move_set.zmax_rot_sum), float(move_set.zmax_rot_sum)]
        
        self.rot_ranges["L"] = [-float(move_set.lmax_rot), float(move_set.lmax_rot)]
        self.rot_ranges_sum["L"] = [-float(move_set.lmax_rot_sum), float(move_set.lmax_rot_sum)]
        
        self.rot_ranges_def = deepcopy(self.rot_ranges_sum)
        
        self.set_translation_vec(move_set.max_trans_vector, move_set.max_trans_vector_sum, simparams)
        
        self.set_move_state()
            
    def set_move_state(self):
        """
        """
        if len(self.allowed_transform) == 0: self.set_state("fixed")
        else: self.set_state("movable")
            
    def set_rot_axis(self, *args):
        """
        """
        for arg in args: self.rot_axis.append(arg)
        
    def set_rotations(self, move_set):
        """
        """
        if move_set.xmax_rot != 0.: self.set_rot_axis("X")
        if move_set.ymax_rot != 0.: self.set_rot_axis("Y")
        if move_set.zmax_rot != 0.: self.set_rot_axis("Z")
        
    def set_state(self, state):
        """
        Arguments:
        ----------
            state  : movable or fixed
        Raises:
        ----------
            Component_movesError for unknown transformation state
        """
        if (state == "fixed") or (state == "movable"): pass
        else: raise Component_movesError("Unknown transformation state %s"%(state))
        self.state = state
            
    def set_transformations(self, *args):
        """
        Arguments:
        ---------
            args : rotation or translation
        Raises:
        ----------
            Component_movesError for unknown transformation
        """
        for arg in args:
            if arg == "rotation" or arg == "translation" or arg == "exchange" or arg == "exchangeandsample" or arg == "SimulateDisorder" or\
                    arg == "rotation_all" or arg == "translation_all" or \
                    arg == "rotation_cov" or arg == "rotation_whole" : pass
            else: raise Component_movesError("Unknown transformation %s"%(arg))
            self.allowed_transform.append(arg)
            
    def set_translation_vec(self, maxtrans, maxtrans_sum, simulparams, default=False):
        
        if not default:
            for i in range(0,3):
                if not str(maxtrans[i]).isalpha(): pass
                elif maxtrans[i].upper() == "X" or maxtrans[i].upper() == "NL": maxtrans[i] = float(simulparams.max_trans_vec[i])
        self.trans_ranges = [[-float(maxtrans[0]), float(maxtrans[0])],\
                             [-float(maxtrans[1]), float(maxtrans[1])], \
                             [-float(maxtrans[2]), float(maxtrans[2])]]
        
        if maxtrans_sum:
            for i in range(0,3):
                if not str(maxtrans_sum[i]).isalpha(): pass
                elif maxtrans_sum[i].upper() == "X" or maxtrans_sum[i].upper() == "NL": maxtrans_sum[i] = 100000000000000000000000000000000000000
            self.trans_ranges_sum = [[-float(maxtrans_sum[0]), float(maxtrans_sum[0])],\
                                    [-float(maxtrans_sum[1]), float(maxtrans_sum[1])],\
                                    [-float(maxtrans_sum[2]), float(maxtrans_sum[2])]]
        
            self.trans_ranges_def = deepcopy(self.trans_ranges_sum)


#@TODO: combine into one function!!

    def set_range_values(self, vmax, vmax_sum, axis, simparams):
        #default rotations
        if str(vmax).isalpha() or str(vmax_sum).isalpha():
            if vmax.upper() == "X" or vmax.upper() == "NL": vmax = float(simparams.max_rot_angle)
            elif vmax_sum.upper() == "X" or vmax_sum.upper() == "NL": vmax_sum = 100000000000000000000000000000000000000
        return vmax, vmax_sum
        
        
class Movable():
    
    def __init__(self, chain_name, line):
        self.xmax_rot = 0.0  #max rot angle in single step
        self.ymax_rot = 0.0
        self.zmax_rot = 0.0
        self.lmax_rot = 0.0 #max rot angle in single step for rotation around line
        
        self.xmax_rot_sum = 0.0  #max rot angle sum in whole simulation
        self.ymax_rot_sum = 0.0
        self.zmax_rot_sum = 0.0
        self.lmax_rot_sum = 0.0
        
        self.max_trans_vector     = 0.0
        self.max_trans_vector_sum = 0.0
        ##call simulated volume components with chain letter identifier only;the same as for pyrystruture.chain
        if "_" not in chain_name: self.comp_name = chain_name
        else:
            self.comp_name = chain_name.split("_")[0]
        line                      = self.process_line(line)
        self.state                = ""    #fixed or movable
        self.set_state(line)
        if self.state != "fixed":
            self.set_rotations(line)
            self.set_rotations_sum(line)
            self.set_translations(line)
            self.set_translations_sum(line)
            
    def process_line(self,line):
        li = line.split()
        line = li[:17]
        #if len(li) == 3 and li[2] == "fixed": return li
        if li[2] == "fixed": return li
        if len(line)<17: raise InputError("You provided to little values for MOVE_STATE %s instead of 16"%(len(line)))
        if not line[-1].isdigit():
            if line[-1] == "NL" or line[-1] == "nl": pass
            elif "." not in line[-1]: raise InputError("You provided to little values for MOVE_STATE %s instead \
                                                    of 16. last value in the raw should be digit, not string"%(len(line)))
        return li[:17]
        
    def set_rotations(self, line):
        """
        user defines what is the maximum rotation angle around xyz axes a IN SINGLE STEP for given component
        """
        self.xmax_rot = self.check_rotation_params([line[3]])
        self.ymax_rot = self.check_rotation_params([line[4]])
        self.zmax_rot = self.check_rotation_params([line[5]])
        self.lmax_rot = self.check_rotation_params([line[15]])
                
    def check_rotation_params(self, values):
        
        rot_sums = []
        for val in values:
            if val.upper() == "X": val=10.
            elif val.upper() == "NL": val = 360.0
            elif val.isalpha(): raise InputError("%s value you provided is not a digit or X or NL"%(val))
            elif float(val)> 360: raise InputError("%s value you provided is not in range 0 to 360"%(val))
            elif float(val)< 0: raise InputError("%s value you provided is not in range 0 to 360"%(val))
            else: val = float(val)
        return val
            
    def set_rotations_sum(self,line):
        """
        user defines what is the maximum sum of rotations around xyz axes for given component IN ALL SIMULATION STEPS
        """
        self.xmax_rot_sum = self.check_rotation_params([line[9]])
        self.ymax_rot_sum = self.check_rotation_params([line[10]])
        self.zmax_rot_sum = self.check_rotation_params([line[11]])
        self.lmax_rot_sum = self.check_rotation_params([line[16]])
        
        #self.check_rotation_params([self.xmax_rot_sum, self.ymax_rot_sum, self.zmax_rot_sum, self.lmax_rot_sum])
        
        if float(self.xmax_rot_sum) < float(self.xmax_rot): raise InputError("Rotation angle in single step cannot be larger than max rotation")
        if float(self.ymax_rot_sum) < float(self.ymax_rot): raise InputError("Rotation angle in single step cannot be larger than max rotation")
        if float(self.zmax_rot_sum) < float(self.zmax_rot): raise InputError("Rotation angle in single step cannot be larger than max rotation")
        if float(self.lmax_rot_sum) < float(self.lmax_rot): raise InputError("Rotation angle in single step cannot be larger than max rotation")
        
                
    def set_state(self, line):
        self.state = line[2]
        
    def check_translation_params(self, values):
        for val in values:
            if val.upper() == "X": val = 10.
            elif val.upper() == "NL": val = 1000000000000000000000000000000000000.0
            elif val.isalpha(): raise InputError("%s value you provided is not a digit or X or NL"%(val))
            else: val = float(val)
        return val
        
    def set_translations(self, line):
        """
        user defines what is the maximum translation vector around xyz axes IN SINGLE STEP for given component
        """
        
        xtrans = self.check_translation_params([line[6]])
        ytrans = self.check_translation_params([line[7]])
        ztrans = self.check_translation_params([line[8]])
        self.max_trans_vector = [xtrans, ytrans,ztrans]
        
    def set_translations_sum(self, line):
        """
        user defines what is the maximum sum of translations around xyz axes for given component IN ALL SIMULATION STEPS
        """
        #self.max_trans_vector_sum = [line[12], line[13],line[14]]
        xtrans = self.check_translation_params([line[12]])
        ytrans = self.check_translation_params([line[13]])
        ztrans = self.check_translation_params([line[14]])
        self.max_trans_vector_sum = [xtrans, ytrans, ztrans]
        
        if float(self.max_trans_vector_sum[0]) < float(self.max_trans_vector[0]): raise InputError("Translation vector in single step cannot be larger than max translation")
        if float(self.max_trans_vector_sum[1]) < float(self.max_trans_vector[1]): raise InputError("Translation vector in single step cannot be larger than max translation")
        if float(self.max_trans_vector_sum[2]) < float(self.max_trans_vector[2]): raise InputError("Translation vector in single step cannot be larger than max translation")
        
