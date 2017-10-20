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

from math                                  import pi
#numpy
from numpy                                 import array
#biopython
from Bio.PDB.Atom                          import Atom
from Bio.PDB                               import NeighborSearch
#PyRy
from Modules.Trans.CCP4_reader             import CCP4
from Modules.Trans.SAXS_reader             import SAXSShape
from Modules.Input.vector                  import count_dist
from Modules.Error.Errors                  import Complex_mapError, InputError
#from Modules.Trans.DensityMapStatistic     import Statistic
from Modules.Trans.ComponentRepresentation import Grid, GridCell


"""@TODO:
alphabetical order of methods
check names of methods and variables
remove code redundancy (has_changed, check_if_changed)
check private vs public methods!!
__str__ method
missing docstrings, descriptions etc
"""

#---------------Minimal Cuboid describing given density map --------------------
#-------------------------------------------------------------------------------
class MapCuboid(object):
    """
        class represents minimal cuboid describing a density map
          _ _ _
        /_| _ _/ |
        | |    | |
        | |_ __|_|
        |/_ __ |/
    """
    
    def __init__(self, xcoord, ycoord, zcoord, radius):
    #cuboid coordinates of min and max points
        self.xmax = 0.000
        self.xmin = 0.000
        self.ymax = 0.000
        self.ymin = 0.000
        self.zmax = 0.000
        self.zmin = 0.000
    #--cuboid main diameter
        self.diameter = 0.000
    #--cuboid center--------------
        self.center = []
    #--lists of all x,y,z coordinates
        self.xcoords = xcoord #e.g. [x1, x2, ....., xn]
        self.ycoords = ycoord
        self.zcoords = zcoord
    #---cuboid edges lengths
        self.a = 0.000
        self.b = 0.000
        self.c = 0.000
	self.point_radius = radius
        
    def __str__(self):
        return "%s %s %s %s %s %s %d %s " % \
          (self.xmax, self.xmin, self.ymax, self.ymin, self.zmax, self.zmin,\
            self.diameter, self.center)
        
    def assign_cuboid_values(self):
        """
             assigns values to mapcuboid attributes
        """
        self.set_vertexes()
        self.set_center()
        self.get_max_distance()
        self.set_edges_lengths()
	
    ###################################################	
    def get_cuboid_coords(self, vertex):
        """
            calculates simulation box size for simbox with no map
        """
	
        self.xmax = self.center[0] + (vertex)/2
        self.ymax = self.center[1] + (vertex)/2
        self.zmax = self.center[2] + (vertex)/2
                
        self.xmin = self.center[0] - (vertex)/2
        self.ymin = self.center[1] - (vertex)/2
        self.zmin = self.center[2] - (vertex)/2
    ###################################################
	
	
    def get_max_distance(self):
        """
            calculates cuboid diameter
        """
        self.diameter = count_dist([self.xmax, self.ymin, self.zmax], \
                                   [self.xmin, self.ymax, self.zmin])
        
    def set_center(self): #Cuboid class
        """
            calculates cuboid center
        """
        self.center = [(self.xmax+self.xmin)/2, (self.ymax+self.ymin)/2, (self.zmax+self.zmin)/2]
    
    def set_edges_lengths(self):
        """
            calculates a,b,c lenghts of cuboid
	    
        """
	self.a = count_dist([self.xmax,self.ymax, self.zmax], #x
                            [self.xmin, self.ymax, self.zmax])
		
        self.b = count_dist([self.xmax,self.ymax, self.zmax],
                            [self.xmax, self.ymin, self.zmax])  #y
	
        self.c = count_dist([self.xmax,self.ymin, self.zmax], #z
                            [self.xmax, self.ymin, self.zmin])
	
    def set_vertexes(self):
        """
            calculates maximum and minimum values of x,y,z coordinates
        """
        self.xmax = max(self.xcoords) 
        self.xmin = min(self.xcoords) 
        self.ymax = max(self.ycoords) 
        self.ymin = min(self.ycoords) 
        self.zmax = max(self.zcoords) 
        self.zmin = min(self.zcoords) 
		
	if self.xmax < 0: self.xmax -= self.point_radius
	else: self.xmax += self.point_radius
	
	if self.xmin < 0: self.xmin -= self.point_radius
	else: self.xmin += self.point_radius
	
	if self.ymax < 0: self.ymax -= self.point_radius
	else: self.ymax += self.point_radius
	
	if self.ymin < 0: self.ymin -= self.point_radius
	else: self.ymin += self.point_radius
	
	if self.zmax < 0: self.zmax -= self.point_radius
	else: self.zmax += self.point_radius
	
	if self.zmin < 0: self.zmin -= self.point_radius
	else: self.zmin += self.point_radius
	

#------ SimulBox describes area where components might be-----------------------
#-------rotated/translated during MC simulation---------------------------------
#-------------------------------------------------------------------------------       
class Simulbox(MapCuboid):
    """
        represents simulation box where complex components will be organized
        in a complex structure
    """
    
    def __init__(self, boxsize):
        self.boxsize     = boxsize #how many times simulation box is bigger than a map
        self.simgrid     = None    #list of Simulbox GridPoints
        #self.ns_simgrid  = []      #NeighborSearch list of all grid cells
        
    def __str__(self):
        return "%s %s %s %s %s %s %s %s" % \
          ( self.xmax, self.xmin, self.ymax, self.ymin, self.zmax, self.zmin,\
            self.center, self.diameter)
	
    def build_simulbox(self, mapcuboid, grid_type, radius):
        """
            assigns values to simulbox attributes
        """
        self.center = mapcuboid.center
        self.get_simulbox_coords(mapcuboid)
              
    def get_simulbox_coords(self, mapcuboid):
        """
            calculates simulation box size
        """
        self.xmax = self.center[0] + (self.boxsize*mapcuboid.a)/2
        self.ymax = self.center[1] + (self.boxsize*mapcuboid.b)/2
        self.zmax = self.center[2] + (self.boxsize*mapcuboid.c)/2
                
        self.xmin = self.center[0] - (self.boxsize*mapcuboid.a)/2
        self.ymin = self.center[1] - (self.boxsize*mapcuboid.b)/2
        self.zmin = self.center[2] - (self.boxsize*mapcuboid.c)/2

	
    def check_spacepoints_locations(self, spacepoints):
	"""
	checks whether spacepoints are inside simulation area, if not,
	they are removed from list of point distance restraints
	"""
	for point in spacepoints:
	    if self.__is_outside(point):
	        print "Point %s is outside simulation area! this restraint\
		       was removed from restraint list"%(point)
	        spacepoints.remove(point)
	return spacepoints
		
    def __is_outside(self, point):
        """
            method checks whether point is within entity
        Arguments:
        ---------
            point - has x,y,z coordinates
        Returns:
        --------
            True if points is outside simulbox
	    False if it is inside
        """
	x, y, z = point.coord[0], point.coord[1], point.coord[2]

        if (self.xmin <= x <= self.xmax) and \
                (self.ymin <= y <= self.ymax) and \
                    (self.zmin <= z <= self.zmax): return False
        else: return True
                
        
#----------Complex_map describes the density map of a whole complex-------------
#-------------------------------------------------------------------------------

class Complex_map(object):
    """
        class that stores information about a density map,  
    """
              
    def __init__(self):
        self.mapfile       = None               # file with a density map if exist
	self.saxsfile      = None               # file with a SAXS shape if exist
        self.center        = None               # actual center of given map
        self.mapcuboid     = None               # minimal cuboid describing given density map
        self.simulbox      = None               # simulbox class object
	self.surfacepoints = []                 #mappoints on the surface of density map
	self.density_sum   = 0.                 #to store sum value of all points' densities
	self.map_threshold = 0.0                #minimal density value describing complex shape
        
    #----experimentally determined points with assigned density
        self.mappoints    = []                  # list of MapPoints' objects
        self.mapgrid      = []                  #mappoints defining volume grid
        self.ns_mapgrid   = []                  #NeighborSearch for all mapgrid cells
	self.map_radius    = 0.0                 #map width from cpp4 file or radius in abinitio saxs reconstruction
        
    #----coordinates-------
        self.xcoordinates = []                  #list of all X coordinates of the grid
        self.ycoordinates = []                  #list of all Y coordinates of the grid
        self.zcoordinates = []                  #list of all Z coordinates of the grid
	
	self.kvol_threshold = 0.0               #stores KVOL parameter value from config file; kvol = 1 means that density map should reflect 1*volume of components structures
	self.grid_type = None                   #type of grid defined by the user in config; at the moment only cubic grid is available
	
    
    def __str__(self):
        return "%s %s %s" % \
          ( self.mapfile, self.center, self.simulbox.center )
        
    def add_mappoint(self, point):
        """
            adds mappoint object into self.mappoints
        """
        self.mappoints.append(point)
	
    def assign_mapcuboid(self, xs, ys, zs, radius, grid_type):
        """
           assigns cuboid describing given density map
        xs, ys, zs - list of x, y, z coordinates
        """
	self.mapgrid = Grid(grid_type, radius, 2*radius, -1.)
	self.grid_type = grid_type
        self.mapcuboid = MapCuboid(xs, ys, zs, radius)
        self.mapcuboid.assign_cuboid_values()
        self.build_mapgrid(self.mapcuboid)
	
	self.mapgrid.set_xyz_ranges(self.mapcuboid.xmin, self.mapcuboid.xmax,\
                                    self.mapcuboid.ymin, self.mapcuboid.ymax,\
                                    self.mapcuboid.zmin, self.mapcuboid.zmax)
	self.mapgrid.generate_cubic_grid(self.mappoints, self.map_radius)
        
	self.mapgrid.set_map_threshold(self.map_threshold)
	self.density_sum = self.mapgrid.get_map_density_sum()
		
	if self.mapgrid.mapcells:
	    self.ns_mapgrid = NeighborSearch(self.mapgrid.mapcells)
	else:
	    self.ns_mapgrid = []
        
    def build_mapgrid(self, mapcuboid):
        self.mapgrid.set_xyz_ranges(mapcuboid.xmin, mapcuboid.xmax,\
                                    mapcuboid.ymin, mapcuboid.ymax,\
                                    mapcuboid.zmin, mapcuboid.zmax)
    
    def create_simulbox(self, simbox_size, grid_type, radius):
        """
           assigns simulation box describing simulation area
Arguments:
----------
    simbox_size : simulbox diameter defined by the user
        """
	
        self.simulbox = Simulbox(simbox_size)
        self.simulbox.build_simulbox(self.mapcuboid, grid_type, radius)
        self.simulbox.get_max_distance()

###################################
#        import scipy.spatial
#	from Modules.Trans.ComponentRepresentation import scipyconverter
#
#        cells, mcells = [],[]
#	
#	
#	scitree = scipyconverter(self.simgrid.grid_cells)
#        self.ns_simgrid = scitree
#	
#	if self.simgrid.mapcells:
#	    scitree = scipyconverter(self.simgrid.mapcells)
#            self.ns_mapgrid = scitree
#	else:
#	    self.ns_mapgrid = []
###################################
	
        
#print "Diam", self.mapcuboid.diameter, self.simulbox.diameter
        
    def get_data_from_ccp4_file(self, simboxradius, threshold = False, kmax = False, first_complex = False):
        """
            calls CCP4_reader and retrieves desired densities
            (of wanted threshold or within assigned volume)
        """
	
	ccp4            = CCP4(self.mapfile)
        count, xcoords, ycoords, zcoords = 0, [], [], []
        index           = 0
	
        if not kmax:
	    densities = ccp4.read_volume_fast(threshold)
	else:
	    print "-----", kmax
	    densities = ccp4.read_simulated_volume_fast(kmax, first_complex)
	    threshold = ccp4.get_map_threshold()
	    self.set_kvol_threshold(threshold)
	
	if len(densities.keys()) == 0:
	    raise Complex_mapError("The density map threshols does not contain any map points with assigned density value. Please change the threshold!!")	
		
	
	
	self.map_threshold = threshold
	self.map_radius = ccp4.width/2.
	if self.map_radius <= 0.: raise InputError ("The cell dimensions in CCP$ file are wrong, should be >0")
	
	for point in densities.keys():
	    x,y,z = point[0], point[1], point[2]
	    rounded_density = round(densities[point], 4)
	    
	    mappoint = GridCell([x,y,z], index)
	    mappoint.set_cell_density(rounded_density)
	    ##mappoint.set_cell_radius(mappoint_radius)
	    self.add_mappoint(mappoint) #list of all MapPoint objects
	    #print "SPHERE, ", str(x), ", ", y, ", ", z, ", ", self.map_radius, ","
	    if x not in xcoords: xcoords.append(x)
	    if y not in ycoords: ycoords.append(y)
	    if z not in zcoords: zcoords.append(z)
	    index += 1
	    
	

	#normalize radii according to densities
	#sorted(density_vals, reverse=True)
	density_vals = ccp4.dotvalues
	
	density_vals = sorted(density_vals, reverse=True)
	distrange = [0.1*self.map_radius, self.map_radius]
	distdif = self.map_radius - 0.1*self.map_radius

	curr = -1.0
	different_densities = 0
	for i in xrange(0, len(density_vals) ):
		if curr != density_vals[i]:
			curr = density_vals[i]
			different_densities += 1

	interval = distdif/different_densities
	
	ms = sorted(self.mappoints, key=lambda mappoint: mappoint.density, reverse=True)
	#print "mappoints", len(ms), len(densities.keys())
	prev = ms[0]
	
	act_radius = self.map_radius
	for mp in ms:	    
	    if mp.density != prev.density:
		mp.radius = act_radius - interval
	    else:
		mp.radius = act_radius
	    prev = mp
	    act_radius = mp.radius
	#print "----", ms[0].radius, ms[0].density
	#print ";;;;;", ms[-1].radius, ms[-1].density
	
        return xcoords, ycoords, zcoords
	
	
    def get_data_from_saxs_file(self, simboxradius, saxsradius):
        """
            calls CCP4_reader and retrieves desired densities
            (of wanted threshold or within assigned volume)
        """
	xcoords, ycoords, zcoords = [], [], []
	
        saxs = SAXSShape()
        saxs.get_data_from_dammif_file(self.saxsfile)
        self.map_radius = saxsradius
	if self.map_radius <= 0.: raise InputError ("You must provide SAXSRADIUS in configuration file")
	
	
        for da in saxs.dummy_atoms:
	    #collect all points within a density map
	    x,y,z = da.coord[0], da.coord[1],da.coord[2]
	    mappoint = GridCell([x,y,z])
	    self.add_mappoint(mappoint)
	    if x not in xcoords: xcoords.append(x)
            if y not in ycoords: ycoords.append(y)
            if z not in zcoords: zcoords.append(z)

	return xcoords, ycoords,zcoords
    
    def get_surface_accesible_mappoints(self):
	"""
	"""
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	#import scipy.spatial
	#from Modules.Trans.ComponentRepresentation import scipyconverter
	
	#scitree = scipyconverter(self.mappoints)
	
        ns = NeighborSearch(self.mappoints)

        #print "----surfacepoints", len(self.surfacepoints)
	for point in self.mappoints:
	    center = array(point.coord)
	    neigh_nr = ns.search(center, self.map_radius*2.05) #2*point.radius)
	    #neigh_nr = scitree.kdtree.query_ball_point(center, 2.*point.radius, eps=0)
	    point.set_neighb_nr(len(neigh_nr))
#@TODO test value of 5 neighbours
	    if len(neigh_nr) < 5:
		self.surfacepoints.append(point)
		
    
    def is_inside(self, point, entity):
        """
            method checks whether point is within entity
        Arguments:
        ---------
            point - has x,y,z coordinates
            entity - object with assigned xmin, xmax, ...., zmin values
        Returns:
        --------
            True if points is inside the enity
            False if it is outside
        """
	if isinstance(point,Atom): x, y, z = point.coord[0], point.coord[1],\
                                                                  point.coord[2]

        #~ print "garnice: sb_xmin: " + str(entity.xmin)
        #~ print "atom: [" + str(x) + ", " + str(y) + ", " + str(z) + "]"
	
        if (entity.xmin <= x <= entity.xmax) and \
                (entity.ymin <= y <= entity.ymax) and \
                    (entity.zmin <= z <= entity.zmax): return True
        else: return False
    
    def set_map_obj(self, simbox_size, grid_type, radius):
        """
            sets object attributes: density, center and mapfile!
        """
        self.create_simulbox(simbox_size, grid_type, radius)
        self.center = self.mapcuboid.center
    
    def get_mappoints(self):
        """
            returns density map points
        """
        for point in self.mappoints: yield point

    def set_mapgrid_by_saxs_data(self, simulbox_radius, saxs_radius, grid_type):
	"""
	    calls methods able to create a grid for piece of
	    dummy atom shape from SAXS experiment
	"""
	xcoords, ycoords, zcoords = self.get_data_from_saxs_file(simulbox_radius, saxs_radius)
	self.assign_mapcuboid(xcoords, ycoords, zcoords, simulbox_radius, grid_type)
	
	######################
	self.get_surface_accesible_mappoints()
	#######################

    def set_mapgrid_by_threshold(self, simulboxradius, threshold, grid_type):
        """
            calls methods able to create a grid for piece of density map
            (depending on user request can be piece of volume or points
            with accurate density value)
Arguments:
---------
    volume : complex volume
    kmax   : how many complex volumes will the density map represent
        """
        xcoords, ycoords,zcoords = self.get_data_from_ccp4_file(simulboxradius, threshold)
        self.assign_mapcuboid(xcoords, ycoords, zcoords, simulboxradius, grid_type)
	
	self.get_surface_accesible_mappoints()

    def set_mapgrid_by_volume(self, simulboxradius, first_cx, kmax, grid_type):
        """
            calls methods able to create a grid for piece of density map
            (depending on user request can be piece of volume or points
            with accurate density value)
Arguments:
---------
    threshold : user defined density map threshold
    (only points with at least this value will be used to describe complex shape)
        """
	
        xcoords, ycoords,zcoords = self.get_data_from_ccp4_file(simulboxradius, kmax = kmax, first_complex=first_cx)
	self.assign_mapcuboid(xcoords, ycoords, zcoords, simulboxradius, grid_type)
		
	######################
	self.get_surface_accesible_mappoints()
	#######################
	
    def set_kvol_threshold(self, threshold):
	self.kvol_threshold = threshold
	
    def set_simbox_with_no_shape_descriptor(self, fcomplex, start_orient, grid_radius, grid_type):
        """
           when modeling is performed with no shape descriptor,
           simulation box is calculated based on complex sequence volume
        """
        ###############sprawdzic jednostke objetosci!!!!!!!!!!!!!!!!!!!!!!!
        vertex = round(fcomplex.volume**(1.0/3),4)
        self.mapcuboid = MapCuboid([],[],[], grid_radius)
        
        
        #if start_orient == True:
            #assign minimal cuboid on current complex
        xcoords, ycoords, zcoords = [], [], []
        for comp in fcomplex.components:
            for atom in comp.pyrystruct.struct.get_atoms():
                if atom.coord[0] not in xcoords: xcoords.append(atom.coord[0])
                if atom.coord[1] not in ycoords: ycoords.append(atom.coord[1])
                if atom.coord[2] not in zcoords: zcoords.append(atom.coord[2])
        self.mapcuboid = MapCuboid(xcoords,ycoords,zcoords, grid_radius)
        self.mapcuboid.assign_cuboid_values()

        self.mapgrid = Grid(grid_type, grid_radius, 2*grid_radius, -1.)
        
        
        self.mapgrid.set_xyz_ranges(self.mapcuboid.xmin, self.mapcuboid.xmax,\
                                    self.mapcuboid.ymin, self.mapcuboid.ymax,\
                                    self.mapcuboid.zmin, self.mapcuboid.zmax)
        
        #self.mapgrid.generate_cubic_grid(self.mappoints, self.map_radius) 
        
        #self.mapgrid.set_map_threshold(self.map_threshold)
        #self.density_sum = self.mapgrid.get_map_density_sum()
                
        if self.mapgrid.mapcells:
            self.ns_mapgrid = NeighborSearch(self.mapgrid.mapcells)
        else:
            self.ns_mapgrid = []
        
