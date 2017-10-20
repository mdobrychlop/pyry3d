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

from numpy                                             import array
from Modules.Error.Errors                              import InteractionError, RestraintError
from Modules.Constans.Constans                         import AMINOACIDS, NUCLEOTIDES
from Modules.Constans.Logfile                          import logfile
#Filtrest3D
from External_Applications.filtrest3D.RestraintsParser import *
from External_Applications.filtrest3D.RestraintsParser import RestraintsParser
from External_Applications.filtrest3D.PModel           import PModel, EmptyStructure

import hashlib
import random

class Pyry_restraint(object):
    
    def __init__(self, res1_name = "", res1_num = None, res2_name = "",\
                       res2_num = None, chain_id = "", atom_name = ""):
        
        self.res1_num  = int(res1_num)   
        self.res1_name = res1_name
        self.res2_num  = int(res2_num)
        self.res2_name = res2_name
        if self.res2_num < self.res1_num: 
            (self.res1_num, self.res2_num) = (self.res2_num, self.res1_num)
        self.chain_id  = chain_id
        self.atom_name = atom_name
        self.restr_id  = ""
        self.score = 0.0
    
    def __str__(self):
        return "%s %s %s %s %s %s"% (self.chain_id, self.atom_name, \
                   self.res1_num, self.res2_num, self.res1_name, self.res2_name)
        
class CommonInteraction(object):
    
    def __init__(self):
        self.uid  = hashlib.md5(str(random.random())).hexdigest()
    
    def assign_components_to_interactions(self, res, components):
        """
            assigns complex_components strucutres into self.components
    Parameters:
    -----------
    res              : restraint object
    components : list of all complex_components objts
    Returns:
    --------
    assigns self.components attribute for a particular restraint
        """
        found_chain = False #to check if chain name from restrfile exists in structures
        for component in components:
            if component.pyrystruct.chain == res.chain_id:
                #check if chain name exist
                if res.chain_id == None: raise InteractionsError("No chain name \
                                        given for residue %s"%(res.res1_name))
                self.components.append(component.pyrystruct)
                self.chains.append(component.pyrystruct.chain)
                #if name of atom is not specified assign CA or C4'
                if res.atom_name == "": self.__assign_default_atoms_types(component, res)
                found_chain = True
                break
        if found_chain == False: raise InteractionError("There is no structure with chain name %s"%(res.chain_id))
        
    def __assign_default_atoms_types(self, component, res):
        """
            assigns default atom name if none is given by the user
        accordingly:
        for proteins CA
        for RNA C4'
        for DNA C4'
        Parameters:
        -----------
            atom_name   : name of particular atom
            mode        : first/second - indicates order of atoms in restraint definition
        """
        if   component.pyrystruct.moltype == "protein"     : res.atom_name = "CA"
        elif component.pyrystruct.moltype == "RNA" or "DNA": res.atom_name = "C4'"
        
        
        
        
        
        
    def check_restraints_correctness(self):
        """
            checks if atoms and residues given in restraints file
            exist in the structure
        Returns:
        --------
            True if atoms exist in structure    
        Raises:
        -------
            RestraintError if the residues or atoms do not exist
        """
        #check if chains exist
        self.find_chain(0, self.first.chain_id)
        self.find_chain(1, self.second.chain_id)
        #check if atoms/residues/res_num are in real structures..
        for i in range(self.first.res1_num, self.first.res2_num+1):
            self.find_entities(0, self.first.chain_id, self.first.res1_name, \
                                                        i, self.first.atom_name)
        for j in range(self.second.res1_num, self.second.res2_num+1):    
            self.find_entities(1, self.second.chain_id, self.second.res1_name, \
                                                       j, self.second.atom_name)
        #check if a particular restraint is not defined within the same structure
        restr_type = self.check_restraints_from_thesame_struc()
        del self.components
        return restr_type
    
    def check_restraints_from_thesame_struc(self):
        """
            checks if restraint is not assigned within one structure
        Raises:
        --------
            TransError if a user provided restraint within one structure
        """
        if self.first.chain_id == self.second.chain_id:
            logfile.write_message("Restraint is given within one structure, it will not be added to score calculation")
            return False
        else: return True #not in the same structure
        
    def find_chain(self, comp_id, chain):
        model =  self.components[comp_id].struct[0]    
        if not model.has_id(chain): raise RestraintError("The structure does not have %s chain"%(chain))
        
    def find_entities(self, comp_id, chain, resname, resnum, atom):
        """
            terminates the program if chain, residue or atom
            provided by the user in the restraints file do not exist in
            corresponding structure  
        Parameters:
        -----------
            comp_id         : component unique identifier
            chain           : chain name
            resname         : residue name
            resnum          : residue number
            atom            : atom name  
        Raises:
        -------
            RestraintError if chain/residue or atom is wrong in given restraint
        """
        model =  self.components[comp_id].struct[0]
        if not model[chain].has_id((' ', resnum, ' ')) or model[chain].has_id(('H_', resnum, ' ')):
            raise RestraintError("The chain %s does not have %s residue"%(chain,resnum))
        if not model[chain][(' ', resnum, ' ')]:
            resname = self.__get_3letter_name(comp_id, resname)
            raise RestraintError("The residue %s does not have %s number.\
                                 Please check %s chain. The other possibility is that sequences are different in FASTA and PDB files.\
                                 Use IDENTIFY_DISORDERS optinion which allows to determine conformations of flexible/disordered regions.\
                                 "%(resname,resnum, chain))
        if atom:
            if not model[chain][(' ', resnum, ' ')].has_id(atom):
                raise RestraintError("The residue%s does not have %s atom. \
    Please check your pdb file or whether representation you have chosen contain this particular atom"%(resnum,atom))
        
    def get_restraints(self, restraint, struc_components):
        """
            extracts data from filtrest3d restraints
        Parameters:
        -----------
            restraint        : single restraint
            struc_components : a list containing all complex_components
        """
                
        for res_id in restraint.reslist1:
            
            if isinstance(res_id, PDBrange): #for range restraints,e.g. (67-73)-(98-105, 108-115)
                self.first = Pyry_restraint(res_id.res1_name,res_id.res1_num, \
                res_id.res2_name,res_id.res2_num,res_id.chain_id,res_id.atom_name)

            elif isinstance(res_id, PDBId):
                self.first = Pyry_restraint(res_id.res_name,res_id.res_num, \
                res_id.res_name,res_id.res_num,res_id.chain_id,res_id.atom_name)
            
            self.assign_components_to_interactions(self.first, struc_components)
            
        for res_id2 in restraint.reslist2:

            if isinstance(res_id2, PDBrange): #for range restraints, e.g. (67-73)-(98-105, 108-115)
                self.second = Pyry_restraint(res_id2.res1_name,res_id2.res1_num, \
                res_id2.res2_name,res_id2.res2_num,res_id2.chain_id,res_id2.atom_name)

            elif isinstance(res_id2, PDBId):
                self.second = Pyry_restraint(res_id2.res_name,res_id2.res_num, \
                res_id2.res_name,res_id2.res_num,res_id2.chain_id,res_id2.atom_name)
            
            self.assign_components_to_interactions(self.second, struc_components)
                
    def __get_3letter_name(self, molid, resname):
        """
            returns 3letter code of given residue  
    Parameters:
    -----------
    molid       : molecule id protein, RNA, DNA
    resname     : residue name
        """
        if self.components[molid].moltype == "protein": return AMINOACIDS[resname]
        else: return NUCLEOTIDES[resname]
        
    def set_di_interaction(self, restraint):
        """
            sets complex_component attributes
        """
        self.restraint = restraint
    
    def set_score(self, score):
        self.score = score
        
    def get_score(self):
        return self.score
    
class SingleComponentInteraction(CommonInteraction):
    
    def __init__(self):
        self.type = "dist"
        super(SingleComponentInteraction, self).__init__()
    
    def get_restraints(self, restraint, struc_components):
        """
            extracts data from filtrest3d restraints
        Parameters:
        -----------
            restraint        : single restraint
            struc_components : a list containing all complex_components
        """
        #self.set_interaction(restraint, surfacepoints)         
                
        for res_id in restraint.reslist1:
            
            if isinstance(res_id, PDBrange): #for range restraints,e.g. (67-73)-(98-105, 108-115)
                self.first = Pyry_restraint(res_id.res1_name,res_id.res1_num, \
                res_id.res2_name,res_id.res2_num,res_id.chain_id,res_id.atom_name)

            elif isinstance(res_id, PDBId):
                self.first = Pyry_restraint(res_id.res_name,res_id.res_num, \
                res_id.res_name,res_id.res_num,res_id.chain_id,res_id.atom_name)
            
            self.assign_components_to_interactions(self.first, struc_components)
            
    def check_restraints_correctness(self):
        """
            checks if atoms and residues given in restraints file
            exist in the structure
        Returns:
        --------
            True if atoms exist in structure    
        Raises:
        -------
            RestraintError if the residues or atoms do not exist
        """
        #check if chains exist
        self.find_chain(0, self.first.chain_id)
        #check if atoms/residues/res_num are in real structures..
        for i in range(self.first.res1_num, self.first.res2_num+1):
            self.find_entities(0, self.first.chain_id, self.first.res1_name, \
                                                        i, self.first.atom_name)
        del self.components
        return 1
    
class DoubleComponentInteraction(CommonInteraction):
    def __init__(self):
        super(DoubleComponentInteraction, self).__init__()
    pass
    
class ResidueDistanceInteraction(DoubleComponentInteraction):
    """
        keeps information about interaction between two complex components
        including restraint distances between atoms or residues within complex components
    """
    def __init__(self, id):
        self.components = []     #list of interacting complex components e.g. [A,B], is deleted after checking restraints
        self.type       = 'dist'     #e.g. 'dist', restraint_type
        self.restraint  = None   #info about weight, dist_value, relation, type, id
        self.first      = None   #first restraint chain
        self.second     = None   #second restraint chain
        self.chains     = []     #contains names of components chains only
        self.id         = id
        super(ResidueDistanceInteraction, self).__init__()
        

########################################################3        
class SymmetryDistances(DoubleComponentInteraction):
    """
    """
    
    def __init__(self, index):
        self.components = []     #list of interacting complex components e.g. [A,B], is deleted after checking restraints
        self.type       = 'symmetry'     #e.g. 'dist', restraint_type
        self.restraint  = None   #info about weight, dist_value, relation, type, id
        self.first      = None   #first restraint chain
        self.second     = None   #second restraint chain
        self.chains     = []     #contains names of components chains only
        self.id         = id
        super(SymmetryDistances, self).__init__()
    
    #def add_symdistance(self, restraint):
    #    self.symmetry_dist.append(restraint)
        
        
class RelationDistances(DoubleComponentInteraction):
    """
    """
    
    def __init__(self, index):
        self.components = []     #list of interacting complex components e.g. [A,B], is deleted after checking restraints
        self.type       = 'relation'     #e.g. 'dist', restraint_type
        self.restraint  = None   #info about weight, dist_value, relation, type, id
        self.first      = None   #first restraint chain
        self.second     = None   #second restraint chain
        self.third      = None
        self.fourth     = None
        self.chains     = []     #contains names of components chains only
        self.id         = id
        super(RelationDistances, self).__init__()
        
    def check_restraints_correctness(self):
        """
            checks if atoms and residues given in restraints file
            exist in the structure
        Returns:
        --------
            True if atoms exist in structure    
        Raises:
        -------
            RestraintError if the residues or atoms do not exist
        """
        #check if chains exist
        self.find_chain(0, self.first.chain_id)
        self.find_chain(1, self.second.chain_id)
        #check if atoms/residues/res_num are in real structures..
        for i in range(self.first.res1_num, self.first.res2_num+1):
            self.find_entities(0, self.first.chain_id, self.first.res1_name, \
                                                        i, self.first.atom_name)
        for j in range(self.second.res1_num, self.second.res2_num+1):    
            self.find_entities(1, self.second.chain_id, self.second.res1_name, \
                                                       j, self.second.atom_name)
        #check if a particular restraint is not defined within the same structure
        restr_type = self.check_restraints_from_thesame_struc()
        
        self.find_chain(2, self.third.chain_id)
        self.find_chain(3, self.fourth.chain_id)
        #check if atoms/residues/res_num are in real structures..
        for k in range(self.third.res1_num, self.third.res2_num+1):
            self.find_entities(2, self.third.chain_id, self.third.res1_name, \
                                                        k, self.third.atom_name)
        for l in range(self.fourth.res1_num, self.fourth.res2_num+1):    
            self.find_entities(3, self.fourth.chain_id, self.fourth.res1_name, \
                                                       l, self.fourth.atom_name)
        #check if a particular restraint is not defined within the same structure
        restr_type = self.check_restraints_from_thesame_struc()
        del self.components
        
        return restr_type      
        
        
    def get_restraints(self, restraint, struc_components):
        """
            extracts data from filtrest3d restraints
        Parameters:
        -----------
            restraint        : single restraint
            struc_components : a list containing all complex_components
        """
                
        for res_id in restraint.reslist1:
            
            if isinstance(res_id, PDBrange): #for range restraints,e.g. (67-73)-(98-105, 108-115)
                self.first = Pyry_restraint(res_id.res1_name,res_id.res1_num, \
                res_id.res2_name,res_id.res2_num,res_id.chain_id,res_id.atom_name)

            elif isinstance(res_id, PDBId):
                self.first = Pyry_restraint(res_id.res_name,res_id.res_num, \
                res_id.res_name,res_id.res_num,res_id.chain_id,res_id.atom_name)
            
            self.assign_components_to_interactions(self.first, struc_components)
            
        for res_id2 in restraint.reslist2:

            if isinstance(res_id2, PDBrange): #for range restraints, e.g. (67-73)-(98-105, 108-115)
                self.second = Pyry_restraint(res_id2.res1_name,res_id2.res1_num, \
                res_id2.res2_name,res_id2.res2_num,res_id2.chain_id,res_id2.atom_name)

            elif isinstance(res_id2, PDBId):
                self.second = Pyry_restraint(res_id2.res_name,res_id2.res_num, \
                res_id2.res_name,res_id2.res_num,res_id2.chain_id,res_id2.atom_name)
            
            self.assign_components_to_interactions(self.second, struc_components)
            
        for res_id3 in restraint.reslist3:
            
            if isinstance(res_id3, PDBrange): #for range restraints,e.g. (67-73)-(98-105, 108-115)
                self.third = Pyry_restraint(res_id3.res1_name,res_id3.res1_num, \
                res_id3.res2_name,res_id3.res2_num,res_id3.chain_id,res_id3.atom_name)

            elif isinstance(res_id3, PDBId):
                self.third = Pyry_restraint(res_id3.res_name,res_id3.res_num, \
                res_id3.res_name,res_id3.res_num,res_id3.chain_id,res_id3.atom_name)
            
            self.assign_components_to_interactions(self.third, struc_components)
            
        for res_id4 in restraint.reslist4:

            if isinstance(res_id4, PDBrange): #for range restraints, e.g. (67-73)-(98-105, 108-115)
                self.fourth = Pyry_restraint(res_id4.res1_name,res_id4.res1_num, \
                res_id4.res2_name,res_id4.res2_num,res_id4.chain_id,res_id4.atom_name)

            elif isinstance(res_id4, PDBId):
                self.fourth = Pyry_restraint(res_id4.res_name,res_id4.res_num, \
                res_id4.res_name,res_id4.res_num,res_id4.chain_id,res_id4.atom_name)
            
            self.assign_components_to_interactions(self.fourth, struc_components)
                
class SurfaceAccessInteraction(SingleComponentInteraction):
    """
    e.g Leu13 "A" is close to density map surface
    """
    def __init__(self, id):
        self.components = []     #list of interacting complex components e.g. [A,B], is deleted after checking restraints
        self.type       = 'access'     #e.g. 'dist', restraint_type
        self.restraint  = None
        self.first      = None   #first restraint chain
        self.chains     = []     #contains names of components chains only
        self.points     = []     #map surface points
        self.id         = id
        super(SurfaceAccessInteraction, self).__init__()

    def __str__(self):               
        return "%s %s %s" %  (self.type, self.first, self.chains)
                            
    def set_sa_interaction(self, restraint, surfacepoints):
        self.restraint = restraint
        self.points    = surfacepoints
    
    
        
class PointDistanceInteraction(SingleComponentInteraction):
    """
    e.g. (10,20,30) Leu13    "D"  < 10A
          point     residue chain distance
    """
    def __init__(self, id):
        self.components = []     #list of interacting complex components e.g. [A,B], is deleted after checking restraints
        self.type       = ''     #e.g. 'dist', restraint_type
        self.restraint  = None
        self.first      = None   #first restraint chain
        self.chains     = []     #contains names of components chains only
        self.points     = []     #X,Y,Z coordinates of point in 3D space
        self.id         = id
        super(PointDistanceInteraction, self).__init__()


    def __str__(self):               
        return "%s %s %s" % ( self.type, self.first, self.chains )
            
    def set_pd_interaction(self, restraint, density_map):
        self.restraint = restraint
        self.points.append(Point(restraint.point))
        self.points    = density_map.simulbox.check_spacepoints_locations(self.points)
        
        
class LogicalRestraint(object):
    def __init__(self, elements):
        self.elements = elements
        chains = []
        for element in self.elements:
            chains += element.chains
        self.chains = list(set(chains))
        self.scores = []
        self.uid = hashlib.md5(str(random.random())).hexdigest()
    
    def update_score(self, scorerFactory):
        self.scores = []
        for element in self.elements:
            scorer = scorerFactory.createScorer(element)
            scorer.score()
            self.scores.append(element.get_score())

class AndRestraint(LogicalRestraint):
    def get_score(self):
        return max(self.scores)
        
class OrRestraint(LogicalRestraint):
    def get_score(self):
        return min(self.scores)
        
class Point(object):
    def __init__(self, point):
        self.coord = array([point[0], point[1], point[2]])
        
class AccessInteraction(ResidueDistanceInteraction):
    """
    e.g Leu13 "A" is close to density map surface
    """
    def __init__(self, id):
        self.components = []     #list of interacting complex components e.g. [A,B], is deleted after checking restraints
        self.type       = ''     #e.g. 'dist', restraint_type
        self.restraint  = None
        self.resi      = None   #first restraint chain
        self.chains     = []     #contains names of components chains only
        self.id         = id

    def __str__(self):               
        return "%s %s %s" %  (self.type, self.first, self.chains)
        
    def check_restraints_correctness(self):
        """
            checks if atoms and residues given in restraints file
            exist in the structure
        Returns:
        --------
            True if atoms exist in structure    
        Raises:
        -------
            RestraintError if the residues or atoms do not exist
        """
        #check if chains exist
        self.find_chain(0, self.resi.chain_id)
        #check if atoms/residues/res_num are in real structures..
        for i in range(self.resi.res1_num, self.resi.res2_num+1):
            self.find_entities(0, self.resi.chain_id, self.resi.res1_name, \
                                                        i, self.resi.atom_name)
        del self.components
        return 1
        
    def get_restraints(self, restraint, struc_components):
        """
            extracts data from filtrest3d restraints
        Parameters:
        -----------
            restraint        : single restraint
            struc_components : a list containing all complex_components
        """
        #self.set_interaction(restraint, surfacepoints)
        
        #print dir(restraint), "\n...", restraint.res.res_num, restraint.id, restraint.highbound, restraint.lowbound, restraint.name, restraint.type        
        
        self.resi = Pyry_restraint(restraint.res.res_name,restraint.res.res_num, \
                restraint.res.res_name,restraint.res.res_num, restraint.res.chain_id,restraint.res.atom_name)

  
        self.assign_components_to_interactions(self.resi, struc_components)
            
    def set_access_interaction(self, restraint):
        self.restraint = restraint
        self.type      = "access"
        #print "::::::", self.restraint.lowbound, self.restraint.highbound
        if float(self.restraint.lowbound) > float(self.restraint.highbound): raise RestraintError("First access value must be smaller than second")
        if float(self.restraint.lowbound) < 0 or float(self.restraint.lowbound) > 100: raise RestraintError("Access values must be in range 0 to 100")
        if float(self.restraint.highbound) < 0 or float(self.restraint.highbound) > 100: raise RestraintError("Access values must be in range 0 to 100")
        
            
