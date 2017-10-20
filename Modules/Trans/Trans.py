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

from Modules.Trans.Component         import Component
from Modules.Trans.Interaction       import ResidueDistanceInteraction, SurfaceAccessInteraction, \
                                            PointDistanceInteraction, SymmetryDistances, LogicalRestraint, \
                                            RelationDistances, AndRestraint, OrRestraint
from Modules.Trans.Complex_map       import Complex_map
from Modules.Error.Errors            import TransError, ConfigError
from Modules.Constans.Logfile        import logfile


class ComplexInteractions(object):
    
    def __init__(self):
        self.sd_interactions = []       #structure distance restraint between 2 Structure Entities
        self.sa_interactions = []       #surface accessibility restraint  some Structure Entities and mappoints from the surface
        self.pd_interactions = []       #point distance restraint between Structure Entity and point in 3D space P(x,y,z)
        self.symmetry_interactions = [] #distances of equal values
        self.relation_interactions = [] #to store restraints describing relation between distances (>, <, =)
        self.logical_interactions = []  #boolean
        self.all = []
    
    def add(self, interaction):
        self.all.append(interaction)
        if isinstance(interaction, SymmetryDistances):
            self.symmetry_interactions.append(interaction)
        elif isinstance(interaction, ResidueDistanceInteraction):
            self.sd_interactions.append(interaction)
        elif isinstance(interaction, SurfaceAccessInteraction):
            self.sa_interactions.append(interaction)
        elif isinstance(interaction, PointDistanceInteraction):
            self.pd_interactions.append(interaction)
        elif isinstance(interaction, RelationDistances):
            self.relation_interactions.append(interaction)
        elif isinstance(interaction, LogicalRestraint):
            self.logical_interactions.append(interaction)
        
        
    def add_sd_interaction(self, interaction):
        self.sd_interactions.append(interaction)
        self.all.append(interaction)
        
    def add_sa_interaction(self, interaction):
        self.sa_interactions.append(interaction)
        self.all.append(interaction)

    def add_pd_interaction(self, interaction):
        self.pd_interactions.append(interaction)
        self.all.append(interaction)
        
    def add_relation_interaction(self, interaction):
        self.relation_interactions.append(interaction)
        self.all.append(interaction)
        
    def add_symmetry_interaction(self, interaction):
        self.symmetry_interactions.append(interaction)
        self.all.append(interaction)

    def add_logical_interaction(self, interaction):
        self.logical_interactions.append(interaction)
        self.all.append(interaction)

class Trans(object):
    """
    class that represents data conversion from data delivered by the user to PyRy objects
    """    
    
    def __init__(self):
        self.components       = []      #keeps all complex components objects
        self.restraints       = None    #keeps restraints data in filtrest format
        self.interactions     = None    #keeps inf about dist_interactions, surface accessibility and point_distance
        self.mapcomponent     = None    #keeps density map as component objt
        
    def __str__(self):
        return "%s %s %s" % (self.components, self.restraints, self.interactions)
        
    def __create_interaction(self, restraint, components, density_map, iteration):
        """
        """
        if restraint.__class__.__name__ == "DistanceRestraint":
            return self.__assign_distant_restraint(restraint, components, iteration)
                
        elif restraint.__class__.__name__ == "SurfaceAccess":
            return self.__assign_surf_access_restraint(restraint, components, density_map.surfacepoints, iteration)

        elif restraint.__class__.__name__ == "PointDistance":
            return self.__assign_pointdist_restraint(restraint, components, density_map, iteration)
            
        elif restraint.__class__.__name__ == "Symmetry":
            return self.__assign_symmetry_restraint(restraint, components, iteration)
        
        elif restraint.__class__.__name__ == "Relation":
            return self.__assign_relation_restraint(restraint, components, iteration)
            
        elif restraint.__class__.__name__ == "AndRestraint":
            return self.__assign_and_restraint(restraint, components, density_map, iteration)
                
        elif restraint.__class__.__name__ == "OrRestraint":
            return self.__assign_or_restraint(restraint, components, density_map, iteration) #?
    
    def __assign_distant_restraint(self, restraint, components, iteration):
        """
        """
        interaction = ResidueDistanceInteraction(iteration)
        interaction.set_di_interaction(restraint)
        interaction.get_restraints(restraint, components)
        if interaction.check_restraints_correctness():
            return interaction
        return None
            
    def __assign_pointdist_restraint(self, restraint, components, density_map, iteration):
        """
        """
        pd_interaction = PointDistanceInteraction(iteration)
        pd_interaction.set_pd_interaction(restraint, density_map)
        pd_interaction.get_restraints    (restraint, components)
        if pd_interaction.check_restraints_correctness():
            return pd_interaction
        return None
                
    def __assign_surf_access_restraint(self, restraint, components, surfacepoints, iteration):
        """
        """
        sa_interaction = SurfaceAccessInteraction(iteration)        
        sa_interaction.set_sa_interaction(restraint, surfacepoints)
        sa_interaction.get_restraints(restraint, components)
        if sa_interaction.check_restraints_correctness():
            return sa_interaction
        return None
    
    def __assign_relation_restraint(self, restraint, components, iteration):
        """
        """
        #print "!!!!assigning RELATIONS"
        relation_interaction = RelationDistances(iteration)        
        relation_interaction.set_di_interaction(restraint)
        relation_interaction.get_restraints(restraint, components)
        if relation_interaction.check_restraints_correctness():
            return relation_interaction
        return None
            
    def __assign_symmetry_restraint(self, restraint, components, iteration):
        """
        """
        #print "!!!!assigning symmetry"
        interaction = SymmetryDistances(iteration)
        interaction.set_di_interaction(restraint)
        interaction.get_restraints(restraint, components)
        if interaction.check_restraints_correctness():
            return interaction
        return None

    def __assign_and_restraint(self, restraint, components, density_map, iteration):
        elements = []
        for rest in restraint.restr_list:
            elements.append( self.__create_interaction(rest, components, density_map, iteration) )
        interaction = AndRestraint(elements)
        return interaction
        
    def __assign_or_restraint(self, restraint, components, density_map, iteration):
        elements = []
        for rest in restraint.restr_list:
            elements.append( self.__create_interaction(rest, components, density_map, iteration) )
        interaction = OrRestraint(elements)
        return interaction
        
#---get complex components and assign them to Complex_component objts-----------
    def get_component(self, pyrystructure, fasta_seq, simul_params, option = None):
        """
            retrieves all structure objects
        Parameters:
        -----------
            pyrystruct      :   chain name 
            fasta_seq       :   sequence of given structure
            grid_type
        """
        component = Component()
        component.set_component(pyrystructure, fasta_seq, simul_params, option)
        self.components.append(component)
        return component
        
    def get_component_no_struct(self, fastaseq, pyrystruct, simul_params):
        """
            encodes sequences with no structure as components
        Parameters:
        -----------
            chain           :   chain name 
            str_seq         :   sequence of given structure
            structure       :   Bio.PDB structure
        """
        component = Component()
        component.set_component_no_struct(fastaseq, pyrystruct, simul_params)
        ###########################################
        self.components.append(component)
        ##########################################
        return component
    
        
    #-----retrieve user restraints and assign them to Dist_interactions objts-------
    def get_interactions(self, restraints, symrestraints, components, density_map):
        """
        retrieves all restraints and assignes them to Distance_interactions objts
        Parameters:
        -----------
            restraints  :
            components   :
        Raises:
        --------
            TransError if a program deletes given restraint
        """
        self.restraints = restraints
        self.interactions = ComplexInteractions()
        index = 0
        for restraint in restraints:
            self.interactions.add( self.__create_interaction(restraint, components, density_map, index) )
            index += 1
            
        if symrestraints:    
            for srestraint in symrestraints:
                #print "here are symetries!!!!!!!!!!!!!!!!!!!!", srestraint.__class__.__name__
                self.interactions.add( self.__create_interaction(srestraint, components, density_map, index) )
                #symmetry_interactions.add_symdistance(restraint)
                #iter += 1
        
        #self.interactions.add_symmetry_interaction(symmetry_interactions)        
        
#----gets mapfile and assigns attributes into Density_map objt------------------
    def get_map(self, saxsfile, mapfile, simul_params, first_complex):
        """
            sets Density_map object and assigns attributes  
        Parameters:
        -----------
            mapfile         : file with density map
            cx_volume       : volume of entire complex
            simul_params    : user defined simulation parameters
            vollist         : list of all complex' residues MWs
        Returns:
        --------
            mapcomponent    : Complex_map object
        """
        mapcomponent = Complex_map()
        #no shape descriptor
        if saxsfile == None and mapfile == None:
            mapcomponent.set_simbox_with_no_shape_descriptor(first_complex, simul_params.start_orient, simul_params.simboxradius, simul_params.simgridtype)
            simul_params.freespace_penalty = [0.0, 0.0]
            #return mapcomponent
        if saxsfile:
            mapcomponent.saxsfile = saxsfile
            #if 0. in simul_params.density_penalty: raise ConfigError ("To score complex compatibility with SAXS shape please use MAP_FREESPACE parameter in a configuration file")
            if simul_params.is_freespace_defined == False:
                print "WARNING! You have not provided penalty value to score a compatibility with complex shape, \
                default weight has been assigned of 1 1"
                logfile.write_file("WARNING! You have not provided penalty value to score a compatibility with complex shape, \
                default weight has been assigned of 1 1\n")
            if 0. in simul_params.freespace_penalty:
                print "WARNING! Bear in mind that a compatibility with complex shape will not be calculated at least \
                in some parts of the simulations. Check MAP_FREESPACE value in configuration file"
                logfile.write_file("WARNING! Bear in mind that a compatibility with complex shape will not be calculated at least \
                in some parts of the simulations. Check MAP_FREESPACE value in configuration file\n")
            
            mapcomponent.set_mapgrid_by_saxs_data( simul_params.simboxradius, simul_params.saxsradius, simul_params.simgridtype)
        #get data when kmax is defined!
        elif simul_params.threshold == None: #if threshold
            if mapfile:
                mapcomponent.mapfile = mapfile
                mapcomponent.set_mapgrid_by_volume(simul_params.simboxradius, first_complex, simul_params.kvol, simul_params.simgridtype)
        elif simul_params.kvol == None:
            if mapfile:
                if 0. in simul_params.freespace_penalty and 0. in simul_params.density_penalty:
                    print "WARNING! Bear in mind that a compatibility with complex shape will not be calculated at least in some parts of the simulations.\
                Check MAP_FREESPACE and DENSITY values in configuration file"
                    logfile.write_file("WARNING! Bear in mind that a compatibility with complex shape will not be calculated at least in some parts of the simulations.\
                Check MAP_FREESPACE and DENSITY values in configuration file")
                if simul_params.is_freespace_defined == False and simul_params.is_density_defined == False: print "WARNING! You have not provided penalty value to score a compatibility with complex shape, default weights have been assigned of 1 1 to MAP_FREESPACE and DENSITY parametrs. \
                                                      Please check whether these values are defined in a configuration file"
                mapcomponent.mapfile = mapfile
                mapcomponent.set_mapgrid_by_threshold(simul_params.simboxradius, simul_params.threshold, simul_params.simgridtype)
        mapcomponent.set_map_obj(simul_params.simbox, simul_params.simgridtype, simul_params.simboxradius)
        return mapcomponent
    
