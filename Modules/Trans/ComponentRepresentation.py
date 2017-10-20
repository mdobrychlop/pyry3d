#!/usr/bin/env python
# -*- coding: utf-8 -*-

#biopython
from Bio                                import PDB
from Bio.PDB                            import PDBParser, PDBIO, NeighborSearch
from Bio.PDB.Atom                       import Atom
from Bio.PDB.Residue                    import Residue
from Bio.PDB.Chain                      import Chain
from Bio.PDB.Model                      import Model
from Bio.PDB.Structure                  import Structure
#glob, math
from glob                               import glob
from math                               import sqrt, fabs
from copy                               import deepcopy
import os
#numpy
from numpy                              import array, arange
#pyry
from Modules.Error.Errors               import ComponentRepresentationError
from Modules.Input.PyRyStructure        import PyRyStructure
from Modules.Constans.Logfile           import logfile
from Modules.Constans.Constans          import MOLWEIGHTS

from numpy import array
#import scipy.spatial

class ComponentRepresentation(object):
    """
    stands for how a given component will be desribed during simulation
    several options are available:
    * fullatom
    * coarse grain models (currently available in biopython are:
        1) CA-Trace/ C4' for nucleic acids
        2) ENCAD 3-point model (CA, O, Side Chain bead)
        3) MARTINI protein model (BB, Side Chain points [S1 to S4]) )
    * grid
        1) cubic
        2) diamond
    """
    
    def __init__(self, comp, repr_type = None):
        self.component = comp #complex component
        self.repr_type = repr_type

    def __str__(self):
        return "%s %s" % (self.component.chain, self.repr_type)
    
    def create_coarse_repr(self, cg_detail = None):
        """
        Creates coarse grain representation of given complex component
        Parameters:
        ----------
            cg_detail   : defines type of coarse grain e.g. CA, ENAD, MARTINI for proteins
        Returns:
        --------
            cg.cg_struct : reduced structure
        """
        cg = CoarseGrain(self.component, cg_detail)
        return cg
    
    def create_grid(self, grid_type, radius = None, overlap = None):
        """
        Creates grid representation of given complex component  
        Parameters:
        -----------
            grid_type   : defines type of grid e.g cubic, diamond
            radius      : indicates size of single grid cell
            overlap     : indicates how much bigger/smaller is the grid than
                          xmin, xmax, ymin, ymax, zmin, zmax
        Returns:
        --------
            grid        : Grid object
        """
        if radius == None: 
            log.write_message("WARNING! You have not provided grid radius, program chose default 1.0A")
            radius = 1.0
        if overlap == None:
            log.write_message("WARNING! You have not provided grid overlap, program chose default 5A")
            overlap = 5.0
        g = Grid(grid_type, radius, overlap)
        g.create_desired_grid()
        #remove old pymol script if exists!
        if os.path.exists(str(self.component.chain)+"pymol.py") == True:
            os.remove(str(self.component.chain)+"pymol.py")
        g.generate_pymol_grid(str(self.component.chain)+"pymol.py")
        return g
    

class CoarseGrain(ComponentRepresentation):
    """
        for proteins:
        * 1pt*c-alpha
        * 2pt*c-alpha / c-beta
        * 3pt*c-alpha / c-beta / side-chain pseudo-centroid OR side-chain centroid
	* sphere - mass center and radius represents components shape
	* ellipsoid - will be added in the future, mass center and X,Y,Z radii represent components shape
    """
    def __init__(self, component_ps, CG_type = 'CA'):
        self.CG_type = CG_type          #stores info about representation type
        self.cg_struct = None           #strores reduced structure as Bio.PDB
        self.fa_struct = deepcopy(component_ps)   #pyrystruct objt strores original (fullatom) structure as Bio.PDB
	self.radius    = 0.0
	self.molweight = 0.0
        self.__set_representation()

    
    def __str___(self):
        pass
    
    def create_ca_representation(self):
        """
        for proteins: Reduces molecule structure to the CA trace only
        for nucleic : Reduces molecule structure to the C4' trace only
        """
        if   self.fa_struct.moltype == 'protein': atom_type = 'CA'
        elif self.fa_struct.moltype == 'RNA'    : atom_type = "C4'"
        elif self.fa_struct.moltype == 'DNA'    : atom_type = "C4'"
                  
        new_struct = self.create_new_chain(self.fa_struct)
        
        self.build_ca_representation(self.fa_struct, new_struct, atom_type)
        
        #self.save_pdb(new_struct, new_struct.id)
        self.cg_struct = new_struct
        return self.cg_struct
    
    def create_cacb_representation(self):
        """
        for proteins: Reduces molecule structure to the CA trace only
        for nucleic : Reduces molecule structure to the C4' trace only
        """
        if   self.fa_struct.moltype == 'protein': atoms_type = ['CA', 'CB']
        elif self.fa_struct.moltype == 'RNA'    : atoms_type = ["C4'", "P"]
        elif self.fa_struct.moltype == 'DNA'    : atoms_type = ["C4'", "P"]
          
        new_struct = self.create_new_chain(self.fa_struct)
        
        self.build_cacb_representation(self.fa_struct, new_struct, atoms_type)
        
        #self.save_pdb(new_struct, new_struct.id)
        self.cg_struct = new_struct
        return self.cg_struct
    
    def create_3p_representation(self):
#TODO: define third atom for aa definition
        """
        for proteins: Reduces molecule structure to the CA, CB, centroid trace only
        for nucleic : Reduces molecule structure to the C4', P, N1/N9 trace only
        """     
        
        new_struct = self.create_new_chain(self.fa_struct)
        if   self.fa_struct.moltype == 'protein':
            self.build_martini_representation(self.fa_struct, new_struct)
        elif self.fa_struct.moltype == 'RNA' or self.fa_struct.moltype == 'DNA':
            self.build_3pnucleic_representation(self.fa_struct, new_struct)
        #self.save_pdb(new_struct, new_struct.id)
        self.cg_struct = new_struct
        return self.cg_struct
    
    
    def create_sphere_representation(self):
        """
	each chain is here represented by centre of mass only
	"""
	new_struct = Structure('sphrere') 
	my_model = Model(0)
	new_struct.add(my_model)
	
	chain_mass_centres, index = [], 1
	my_chain = Chain(self.fa_struct.chain)
	new_struct[0].add(my_chain)
	    
	coord, self.molmass, self.radius = self.calculate_centre_of_complex(self.fa_struct.struct)
	my_residue = Residue((' ',index,' '),"ALA",' ')
	    
	coords = array(coord,'f')
	atom = Atom('CA',coords, 0, 0, ' ', ' CA',1)
    
	my_chain.add(my_residue)
	my_residue.add(atom)

	self.cg_struct = new_struct
	name = "dddd"+self.fa_struct.chain
	self.save_pdb(new_struct, name)
    
    def calculate_centre_of_complex(self, component):
        """
        calculates centre of mass for the whole complex
        """

        component_centre, max_coordinates = [0.,0.,0.], []

        total_mass = 0

        for atom in component.get_atoms():
	    max_coordinates.extend([fabs(atom.coord[0]), fabs(atom.coord[1]), fabs(atom.coord[2])])
            mass = self.assign_molweight(atom.get_name())
	    total_mass += mass
	    
            component_centre += atom.coord * mass
        radius = max(max_coordinates)
        component_centre /= total_mass

        return component_centre, total_mass, radius
    
    def assign_molweight(self, atom_id):
        """
            assignes a molecular weight to each atom in a given structure
        Raises:
        -------
            Cmplx_ComponentsError if atom name is not known
        """
	for char in atom_id:
	    if char in MOLWEIGHTS.keys():
		atom_name = char
		break
		
        if atom_name in MOLWEIGHTS.keys(): 
            molweight = MOLWEIGHTS[atom_name]
            return molweight
        else: raise PyRyStructureError("Atom not known"+atom_name)

    
    def save_pdb(self, struct, name):
        out = PDBIO()
        out.set_structure(struct)
        out.save(str(name)+'bbbbbb.pdb')
    
#TODO jako argument funkcji podac ktore atomy sa w tej reprezentacji!!!
                
    def build_ca_representation(self, old_struct, new_st, atom_name):
        """
           builds new structure composed of atom_name atoms only
           for nucleic acids these are C4', for proteins CA
        """
        for atom in old_struct.struct.get_atoms():
            if atom.name == atom_name:
                resi = atom.get_parent()
                new_resi = Residue(resi.id, resi.resname, " ")
                new_chain = list(new_st.get_chains())[0]
                new_chain.add(new_resi)
                new_resi.add(atom)
        del self.fa_struct
        
    def build_cacb_representation(self, old_struct, new_st, atom_names):
        """
           builds new structure composed of atom_name atoms only
           for nucleic acids these are ["C4'", "P"], for proteins ["CA", "CB"]
        """
        for atom in old_struct.struct.get_atoms():
            new_chain = list(new_st.get_chains())[0]
            if atom.name in atom_names:
                resi = atom.get_parent()
                if new_chain.has_id(resi.id): pass
                else:
                    new_resi = Residue(resi.id, resi.resname, " ")
                    new_chain.add(new_resi)
                new_resi.add(atom)
        del self.fa_struct
        
    def build_martini_representation(self, old_struct, new_st):
        """
        Reduces complexity of protein residue to the MARTINI coarse grained model:
        CA, O, Bead(s) in specific atom location.
        Reference:
        
        Monticelli et al. The MARTINI coarse-grained force field: extension to proteins.
        J. Chem. Theory Comput. (2008) vol. 4 (5) pp. 819-834
        """
        
        conversion_table = {
                        "L": ["CG"], "V": ["CB"], "E": ["CB"], "d": ["CG"],
                        "T": ["CB"], "Q": ["CB"], "S": ["CB"],  
                        "A": [], "G": [],  
                        "K": ["CG", "CE"],
                        "I": ["CD1"],
                        "R": ["CG", "NE"],
                        "N": ["CG"], "P": ["CG"],
                        "F": ["CG", "CE1", "CE2"], "Y": ["CG", "CE1", "CE2"],
                        "H": ["CB", "ND1", "NE2"],
                        "M": ["CG"],
                        "W": ["CB", "CD1", "CD2", "CE2"],
                        "C": ["SG"],
			    
			"LEU": ["CG"], "VAL": ["CB"], "GLU": ["CB"], "ASP": ["CG"],
                        "THR": ["CB"], "GLN": ["CB"], "SER": ["CB"],  
                        "ALA": [], "GLY": [],  
                        "LYS": ["CG", "CE"],
                        "ILE": ["CD1"],
                        "ARG": ["CG", "NE"],
                        "ASN": ["CG"], "PRO": ["CG"],
                        "PHE": ["CG", "CE1", "CE2"], "TYR": ["CG", "CE1", "CE2"],
                        "HIS": ["CB", "ND1", "NE2"],
                        "MET": ["CG"],
                        "TRP": ["CB", "CD1", "CD2", "CE2"],
                        "CYS": ["SG"],
                        }
		    
        for resi in old_struct.struct.get_residues():
            new_chain = list(new_st.get_chains())[0]
	    #print "@@@", resi.resname, resi.id
	    hetname = resi.id[0].split("_")
	    if hetname[0] == "H": continue
            for atom in resi:
                resname = resi.resname.strip()
		
                if atom.name in conversion_table[resname] or atom.name == "CA":
                    resi = atom.get_parent()
                    if new_chain.has_id(resi.id): pass
                    else:
                        new_resi = Residue(resi.id, resi.resname, " ")
                        new_chain.add(new_resi)
                    new_resi.add(atom)
        del self.fa_struct
        
    def build_3pnucleic_representation(self, old_struct, new_st):
        """
           builds new structure composed of atom_name atoms only
           for nucleic acids these are ["C4'", "P", "N1"/"N9"]
           N1 is taken for ADE and GUA; N9 for CYT, URA, THY
        """
        for resi in old_struct.struct.get_residues():
	    if resi.has_id("N9"):
		atom_names = ["C4'", "P", "N9"]
	    else:
		atom_names = ["C4'", "P", "N1"]
            new_chain = list(new_st.get_chains())[0]
	    for atom in resi:
		if atom.name in atom_names:
		    resi = atom.get_parent()
		    if new_chain.has_id(resi.id): pass
		    else:
			new_resi = Residue(resi.id, resi.resname, " ")
			new_chain.add(new_resi)
		    new_resi.add(atom)
        del self.fa_struct
    
    def create_new_chain(self, old_struct):
        s = Structure(old_struct.chain) 
        my_model = Model(0)
        s.add(my_model)
        my_chain = Chain(old_struct.chain)
        my_model.add(my_chain) #what if more chains in one component?
        return s   
    
    def __set_representation(self):
        """
            
        """
        if self.CG_type == "ca":
            self.create_ca_representation()
        elif self.CG_type == 'cacb':
            self.create_cacb_representation()
        elif self.CG_type == '3p':
            self.create_3p_representation()
	elif self.CG_type == 'sphere':
            self.create_sphere_representation()
	elif self.CG_type == 'ellipsoid':
            self.create_ellipsoid_representation()
        else: raise ComponentRepresentationError('No such coarse grain representation \
                                                        type %s'%(self.CG_type))
    
    
class Grid(ComponentRepresentation):
    """
        abstract class; container for cubicgrid and diamondgrid
    """
    
    def __init__(self,  grid_type, radius, overlap, defscore = None):
        self.grid_type     = grid_type      #cubic, diamond etc.
        self.radius        = radius            #single grid cell size
        self.diameter      = 2*radius*sqrt(3) #single cell diameter
        self.overlap       = overlap          #size of overlap in grid according to structure
        self.grid_cells    = []            #all cells describing given grid
        self.default_score = defscore   #default score value for all grid cells
        self.mapcells      = []              #stores grid cells with map elements inide
	self.mapdensity    = 0.
	self.min_density_value = 1000000000000000          #smallest density occuring in the map
        
    def __str__(self):
        return "%s %s %s %d"%(self.grid_type, self.radius, self.overlap, len(self.grid_cells))
        
    def add_grid_point(self, gridcell):
        if isinstance(gridcell, GridCell): self.grid_cells.append(gridcell)
    
    def calculate_boundaries(self, struct):
        """
        calculates grid boudaries for Bio.PDB structure ranges are defined by:
            xmin, xmax,
            ymin, ymax,
            zmin, zmax
        """
        xcoords, ycoords, zcoords = [], [], []
        for at in struct.get_atoms():
            xcoords.append(at.coord[0])
            ycoords.append(at.coord[1])
            zcoords.append(at.coord[2])
        self.xmax = max(xcoords)
        self.xmin = min(xcoords)
        self.ymax = max(ycoords)
        self.ymin = min(ycoords)
        self.zmax = max(zcoords)
        self.zmin = min(zcoords)
        
    def clear_grid(self):
        """
           assign -1 score to all grid cells
        """
        if len(self.grid_cells) == 0: raise ComponentRepresentationError("Cannot clear empty Grid")
        for point in self.get_cells(): point.score = -1.0
    #
    def generate_cubic_grid(self, mappoints, pixel_size = False, normalization = False):
        """
           builds cubic grid and finds mappoints within it; all gridcells
           containing mappoint have ismap of value True assigned
        """
	#print "=====", pixel_size
	
        if mappoints:
	    ns = NeighborSearch(mappoints) #mappoints must have coord attribute!!!
	    if pixel_size < self.radius: rad = self.diameter/2.0 
	    else: rad = self.diameter/2.0 + pixel_size 
        self.get_xyz_grid_coordinates()
        index = 0                      #number of given gridcell
        for x in self.xcoords:
            for y in self.ycoords:
                for z in self.zcoords:
                    gc = GridCell([x,y,z], index, self.default_score)
                    index += 1
                    #here call NS search through all mappoints radius = simboxradius or half diameter!!
                    if mappoints:
		        center = array([x, y, z])		
			neighbours = ns.search(center, rad, 'A')
			if len(neighbours) != 0:
			    if not normalization:
				gc.ismap = True
			    gc.set_cell_radius(self.radius)
			    #density = gc.set_averaged_density(neighbours)
			    if not normalization:
				self.mapcells.append(gc)			
                    self.grid_cells.append(gc)
	#print "!!!!", len(self.mapcells), len(self.grid_cells)
		    
#	if mappoints:
#	    scitree = scipyconverter(mappoints) #mappoints must have coord attribute!!!
#	    if pixel_size <= self.radius: rad = self.diameter/2.0 
#	    else: rad = self.diameter/2.0 + pixel_size
#	    
#        self.get_xyz_grid_coordinates()
#        index = 0                      #number of given gridcell
#        for x in self.xcoords:
#            for y in self.ycoords:
#                for z in self.zcoords:
#                    gc = GridCell([x,y,z], index, self.default_score)
#                    index += 1
#                    #here call NS search through all mappoints radius = simboxradius or half diameter!!
#                    if mappoints:
#		        center = array([x, y, z])		
#			neighbours = scitree.kdtree.query_ball_point(center, rad, eps=0)
#			if len(neighbours) != 0:
#			    if not normalization:
#				gc.ismap = True
#			    gc.set_cell_radius(self.radius)
#			    #density = gc.set_averaged_density(neighbours)
#			    if not normalization: self.mapcells.append(gc)			
#                    self.grid_cells.append(gc)		    
		    
    #if mappoint radiuses are normalized use normalization = True:
	if normalization and mappoints:
	    ns2 = NeighborSearch(self.grid_cells)	    
	    for mp in mappoints:
		cent = array(mp.coord)
		if mp.radius <= self.radius: r = self.diameter/2.0 #((self.diameter)/2.) + mappoints[0].radius/2.
		else: r = self.diameter/2.0 + mp.radius
		nei = ns2.search(cent, r, 'A')
		if len(nei) != 0:
		    for n in nei:
			n.ismap = True
			n.set_cell_radius(self.radius)
			if n not in self.mapcells: self.mapcells.append(n)
			
	#if normalization and mappoints:
	#    scitree = scipyconverter(self.grid_cells)
	#    for mp in mappoints:
	#	cent = array(mp.coord)
	#	if mp.radius <= self.radius: r = self.diameter/2.0 #((self.diameter)/2.) + mappoints[0].radius/2.
	#	else: r = self.diameter/2.0 + mp.radius
	#	nei = scitree.kdtree.query_ball_point(cent, r, eps=0)
	#	if len(nei) != 0:
	#	    for n in nei:
	#		self.grid_cells[n].ismap = True
	#		self.grid_cells[n].set_cell_radius(self.radius)
	#		if self.grid_cells[n] not in self.mapcells: self.mapcells.append(self.grid_cells[n])			
		    
    def get_map_density_sum(self):
		
	for mappoint in self.mapcells:
	    if mappoint.density == self.min_density_value:
		mappoint.density = 0
	    else:
	        self.__rescale_density_values(mappoint, self.min_density_value)
	    self.mapdensity += mappoint.density
	return self.mapdensity
    
    def set_map_threshold(self, density):
	self.min_density_value = density
    
    def __rescale_density_values(self, mappoint, min_dens):
	"""
	"""
	mappoint.density += -(min_dens)
                    
    def generate_diamond_grid(self):
        pass
        
    def generate_pymol_grid(self, filename):
        """
            prepares python script which is drawing a grid in pymol
        """
        fh = open(filename, "a")
        fh.write("from pymol.cgo import *\nfrom pymol import cmd \n\n")
        fh.write("obj = [\n    BEGIN, LINES, \n    COLOR, 1.0, 1.0, 1.0,\n")
        for cell in self.grid_cells:
            fh.write("    SPHERE, "+str(cell.center[0])+" , "+str(cell.center[1])\
                                +" , "+str(cell.center[2])+" , "+str(1.0)+",\n")
        fh.write("END\n    ] \ncmd.load_cgo(obj,'cgo01')")
        fh.close()
    
    def get_cell_diameter(self):
        return self.diameter
        
    def get_cell_radius(self):
        return self.radius
    
    def get_cells(self):
        """
            returns all cells of given grid
        """
        for cell in self.grid_cells: yield cell
    
    def get_xyz_grid_coordinates(self):
        """
            genrates x, y, z coordinates lists for grid generation process
        """
        self.xcoords = arange(self.xmin - self.overlap, \
                              self.xmax + self.overlap, 2 * self.radius)
        self.ycoords = arange(self.ymin - self.overlap, \
                              self.ymax + self.overlap, 2 * self.radius)
        self.zcoords = arange(self.zmin - self.overlap, \
                              self.zmax + self.overlap, 2 * self.radius)
        
    def set_xyz_ranges(self, xmin, xmax, ymin, ymax, zmin, zmax):
        """
        """
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.zmin = zmin
        self.zmax = zmax
                    
    ##############################################################
#"""    
#    add method to iterate though all grid cells and assigns MAP
#     value to one that contain MAPPOINT!!
#"""    
    #############################################################
                    
class CubicAtomGrid(Grid):

    """
       Grid made if identical cuboids
    """
    def __init__(self, grid):
        self.grid = grid                                  #main "abstract" grid
        self.grid_cells = []                              #stores all grid cells
        self.used_atoms = []                              #helping list
        self.diameter = 2.* self.grid.radius * sqrt(3)    #single cell diameter
#--------x,y,z coordinates for grid genration-------------
   
    def __str__(self):
        return "%s %s"%(len(grid_cells), self.diameter)
    
    
class DiamentGrid(Grid):
    """
       Grid made if identical diamonds
    """

class GridCell(object):
    """
        describes single cell of particular grid; each cell is defined
        by center ([x,y,z] coordinates) and atoms which are inside (PDB objts)
    """
    def __init__(self, center, index=None, score=None):
        self.coord = center            #cell center
        self.content = []              #list of objects inside this cell
        self.score = score             #can indicate score, number of objects, weight etc
        self.density = 0.0
        self.ismap = False             #is given cell a mappoint or not
        self.index = index
        self.radius = 0.
	self.neighb_nr = 0             #number of neighbouring cells/points
        
    def __str__(self):
        return "%s %s %s"%(self.coord, self.radius, self.density)
        
    def add_entity(self, atoms):
        self.content.append(atoms)
    
    def change_score(self, new_score):
        if new_score == 0 or new_score == -1: pass
        else: raise ComponentRepresentationError("Score for grid cel can be 0 or -1!")
        self.score = new_score
	
    def set_averaged_density(self, neighbours):
	densities = 0.
	for neigh in neighbours:
	    densities += neigh.density
	    
	self.density = densities/len(neighbours)
	return self.density
        
    def get_coord(self):
        return self.coord
    
    def set_cell_density(self, density):
        self.density = density
        
    def set_cell_radius(self, rad):
        self.radius = rad
	
    def set_neighb_nr(self, neighb_nr):
	self.neighb_nr = neighb_nr
	
#class scipyconverter():
#    def __init__(self, oblist):
#	self.oblist = oblist
#	self.kdtree = None
#	self.coords = []
#	self.get_coords()
#	self.set_kdtree()
#	#self.ns()
#	
#    def get_coords(self):
#	for a in self.oblist:
#	    self.coords.append(a.coord)
#	self.coords = array(self.coords)
#	
#    def ns(self):
#	for a in self.oblist:
#	    res = self.run_ns(a.coord)
#	print res
#	print res[0].coord, res[0].index
#    
#    def set_kdtree(self):
#	self.kdtree = scipy.spatial.cKDTree(self.coords)
#	
#    def run_neisearch(self, point, r, eps=0):
#	res = []
#        bla = self.kdtree.query_ball_point(point, r)
#	for b in bla:
#	    res.append(self.oblist[b])
#	return res


if __name__=='__main__':

    struct_vol_list = [] #list of all complex components
    #calculate volume for each component
    grid_radius = 3.0
    grid_overlap = 5.0
    for fileh in glob("*.pdb"):
        print "structure", fileh
        p = PDBParser(PERMISSIVE=False, QUIET=True)
        st = p.get_structure('r',fileh)
        gg = CubicGrid(st, grid_radius, grid_overlap)
        gg.generate_grid()
