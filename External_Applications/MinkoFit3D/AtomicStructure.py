__author__ = "Wojtek Potrzebowski"
__copyright__ = "Copyright 2009"
__credits__ = ["Janusz Bujnicki"]
__license__ = "GPL"
__version__ = "0.1.0"
__maintainer__ = "Wojtek Potrzebowski"
__email__ = "wpotrzeb@genesilico.pl"
__status__ = "Prototype"

from Bio.PDB import PDBParser
from Bio.PDB import PDBIO
from Bio.PDB.Structure import Structure
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from numpy import array
from atomic_masses import atomic_mass,  residue_mass,  residue_volume

class AtomicStructure():
    def __init__(self):
        self.pdb_file = None
        self.weighted_coordinates = None
        self.to_file_coordinates = None
        self.protein_psi = None
        self.subunit_pdb =None
        self.coordens = None
        self.molecule_mass = 0.0
        self.atomic_weights = []
        self.pdb_faces = 0
        #Estimated threshold value calulated when
        #strcuture is convereted to mesh
        self.estimated_threshold  = 0.5
        self.p=PDBParser()
    
    def read(self,  subunit_pdb_file):
        """
        Reading structures coordinates - keeping in the init wasn't 
        very good idea
        """
        self.pdb_file = subunit_pdb_file
        self.subunit_pdb = self.p.get_structure("STR", subunit_pdb_file)
        self.__read_subunit_coordinates()
    
    
    #This is specific for definer I will have to rename  it
    def __read_subunit_coordinates(self):
        """
        Reading coordinates from the subunit pdb file
        Both weighted and all coordiantes are read
        """
        coordinates = []
        all_coordinates = []
        density_coordinates = []
        for model in self.subunit_pdb:
            for chain in model:
                for residue in chain:
                    residue_center = [0.0,  0.0,  0.0]
                    total_weight = 0.0
                    for atom in residue:
                        #TODO: Change for nucleic acids here
                        weight = atomic_mass[atom.name[0]]
                        coords = atom.get_coord()
                        all_coordinates.append(coords)
                        #If CA model is read then different weight has to be assaiged for mesh production
                        if len(residue)==1 and atom.name[0] == "C":
                            #Roughly calculated  for glycine resiude - seems to provide reasonable results
                            weight = 37
                        density_coordinates.append([coords[0],coords[1], coords[2], weight])
                        #Calculating geometrical center of the resiude 
                        #In case of CA model the coordinates are not changed
                        residue_center+=coords*weight
                        total_weight +=weight
                        self.atomic_weights.append(weight)
                    residue_center = residue_center/total_weight
                    #This is backpup for nonstandard residudes
                    #The same mass is calculated from the the atomic mass
                    #TODO: Not very neat
                    try:
                        total_weight = residue_mass[residue.resname]
                    except:
                        pass
                    self.molecule_mass+=total_weight
                    #Here we will add center of the residue
                    coordinates.append(residue_center)
        self.weighted_coordinates = array(coordinates)
        self.to_file_coordinates = array(all_coordinates)
        self.coordens = array(density_coordinates)
