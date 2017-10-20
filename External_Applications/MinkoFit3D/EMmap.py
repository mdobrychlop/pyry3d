__author__ = "Wojtek Potrzebowski"
__copyright__ = "Copyright 2009"
__credits__ = ["Janusz Bujnicki"]
__license__ = "GPL"
__version__ = "0.1.0"
__maintainer__ = "Wojtek Potrzebowski"
__email__ = "wpotrzeb@genesilico.pl"
__status__ = "Prototype"

#This class is almost directly taken from situs.
#The appropriate copyright statement will be necessary here
from ccp4_reader import CCP4
from Bio.PDB.NeighborSearch import NeighborSearch
from numpy import array
from numpy import zeros
from numpy import where
from numpy import reshape
from copy import deepcopy

class Neighbor():
    """
    Atom class is required for the niegbor search algorithm.
    It stores coordinates of the atom. 
    """
    def __init__(self,  coordinates):
        """
        Coordinates are given as list
        """
        self.coord=array(coordinates)
    def get_coord(self):
        return self.coord


class EMmap():
    def __init__(self, vol_file,  threshold,  resolution=None):
        self.volume_file=vol_file
        self.threshold = threshold
        self.resolution = resolution
        self.density_coords = {}
        self.simulated_density_coords = {}
        self.density_coords_list = []
        #self.protein_phi = None
        self.map_phi = None
        self.gaussian_phi = None
        self.neighbor_list = []
        self.neighbor_search = None
        self.centroids = []
        self.extx = 0
        self.exty =  0
        self.extz =  0
        self.widthx = 0
        self.widthy = 0
        self.widthz = 0
        self.width = 0
        self.gridx = 0
        self.gridy = 0
        self.gridz = 0
        #Denisties over threshold values
        self.dotvalues = []
        self.map_mean = 0
    
    def read_volume_fast(self):
        """
        Reading volume and returning dictionary [x,y,z]=density
        """
        #CCP4 file reader
        ccp4 = CCP4(self.volume_file)
        densities = ccp4.return_density_values()
        v = ccp4.return_header_values()
        nx, ny, nz = v['nx'], v['ny'], v['nz']
        mx, my, mz = v['mx'], v['my'], v['mz']
        mxstart, mystart, mzstart= v['mxstart'], v['mystart'], v['mzstart']
        xlen, ylen, zlen = v['xlen'], v['ylen'], v['zlen']
        self.extx = nx
        self.exty = ny 
        self.extz = nz 
        self.widthx = xlen/mx
        self.widthy = ylen/my
        self.widthz = zlen/mz
        self.width = max(self.widthx, self.widthy , self.widthz)
        self.gridx = mxstart * self.widthx
        self.gridy = mystart * self.widthy
        self.gridz = mzstart * self.widthz
        #Format independent part starts here
        nvox = self.extx * self.exty * self.extz
        self.map_phi = zeros(nvox)
        #self.protein_phi = zeros(nvox)
        extxy = self.extx * self.exty
        densities = reshape(densities, nvox)
        counts = where(densities>self.threshold)
        for count in counts[0]:
            indv = int(count)
            indz = indv /extxy
            indv = indv - indz * extxy
            indy = indv / self.extx
            indx = indv - indy * self.extx
            x = float(indx)*self.widthx+self.gridx
            y = float(indy)*self.widthy+self.gridy
            z = float(indz)*self.widthz+self.gridz
            #Here there will be the changes for the coordinates and density values
            density_value = densities[count]
            self.density_coords[(float(x), float(y), float(z))] = density_value
            #Needed for neighboor search I guess there should be interaface between EMmap and AtomicStrcucture  class
            self.neighbor_list.append(Neighbor([float(x), float(y), float(z)]))
            self.dotvalues.append(density_value)
            self.map_phi[count]+=float(density_value)
            #TODO: Check this!!!
            ##count = count +1
        self.density_coords_bkp = deepcopy(self.density_coords)
    
    def read_simulated_volume_fast(self,  densities,  structure):
        """
        Reading volume and returning dictionary [x,y,z]=density
        """
        #TODO: In the mesher this values are reread. It can be done once
        ccp4 = CCP4(self.volume_file)
        v = ccp4.return_header_values()
        nx, ny, nz = v['nx'], v['ny'], v['nz']
        mx, my, mz = v['mx'], v['my'], v['mz']
        mxstart, mystart, mzstart= v['mxstart'], v['mystart'], v['mzstart']
        xlen, ylen, zlen = v['xlen'], v['ylen'], v['zlen']
        self.extx = nx
        self.exty = ny 
        self.extz = nz 
        self.widthx = xlen/mx
        self.widthy = ylen/my
        self.widthz = zlen/mz
        #strcuture_dimensions=structure.return_structure_span()
        self.width = max(self.widthx, self.widthy , self.widthz)
        self.gridx = mxstart * self.widthx
        self.gridy = mystart * self.widthy
        self.gridz = mzstart * self.widthz
        #Format independent part starts here
        nvox = self.extx * self.exty * self.extz
        extxy = self.extx * self.exty
        densities = reshape(densities, nvox)
        #############################
        molecule_mass = structure.molecule_mass
        Nv=float(molecule_mass)/(0.81*self.widthx*self.widthy*self.widthz)
        Nv = int(Nv*1.2)
        dot = densities[densities>0.0]
        dot = sorted(dot)
        threshold = round(dot[-Nv], 2)
        structure.estimated_threshold = threshold
        #############################
        print "<<DEBUG: Threshold of simulated map is set to",  threshold
        counts = where(densities>threshold)
        for count in counts[0]:
            indv = int(count)
            indz = indv /extxy
            indv = indv - indz * extxy
            indy = indv / self.extx
            indx = indv - indy * self.extx
            x = float(indx)*self.widthx+self.gridx
            y = float(indy)*self.widthy+self.gridy
            z = float(indz)*self.widthz+self.gridz
            #Here there will be the changes for the coordinates and density values
            density_value = densities[count]
            self.simulated_density_coords[(float(x), float(y), float(z))] = density_value
            self.dotvalues.append(density_value)
            count = count +1

