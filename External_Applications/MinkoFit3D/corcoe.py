__author__ = "Wojtek Potrzebowski"
__copyright__ = "Copyright 2009"
__credits__ = ["Janusz Bujnicki"]
__license__ = "GPL"
__version__ = "0.1.0"
__maintainer__ = "Wojtek Potrzebowski"
__email__ = "wpotrzeb@genesilico.pl"
__status__ = "Prototype"

from Bio.PDB import PDBParser, PDBIO
from numpy import array
from numpy import dot
from numpy import sqrt
from numpy import zeros
from numpy import mean
from numpy.linalg import norm
from ccp4_reader import CCP4
from math import floor, ceil
#from EMmap_sim import EMmap
from AtomicStructure import AtomicStructure
import sys

#This class uses  code  first introduced in Situs.
#Willy Wriggers. Using Situs for the Integration 
#of Multi-Resolution Structures. 
#Biophysical Reviews, 2010, Vol. 2, pp. 21-27

class CorCoe():
    def __init__(self, volume,  structure):
        self.volume=volume
        self.structure = structure
        self.protein_phi = None
        self.map_phi = None
    
    def __generate_protein_map(self):
        """
        Main function for map geneartion
        """
        
        self.protein_phi = zeros((self.extx)*(self.exty)*(self.extz))
        #No mass measure mode - weight = 1 - it might be nice solution too
        weight = 1.0
        #Coordinates will be taken directly from the pdb
        for model in self.pdb_object:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        #Bfactor thresholding
                        #if atom.get_bfactor()<self.b_threshold:
                            coords = atom.get_coord()
                            # compute position within grid - width will be taken from the CCP4 file
                            gx = (coords[0]-self.gridx)/self.widthx
                            gy = (coords[1]-self.gridy)/self.widthy
                            gz = (coords[2]-self.gridz)/self.widthz
                            x0 = floor(gx)
                            y0 = floor(gy)
                            z0 = floor(gz)
                            x1 = x0+1
                            y1 = y0+1
                            z1 = z0+1
                            #interpolate
                            a = x1-gx
                            b = y1-gy
                            c = z1-gz                            
                            reorient = lambda k, j, i : int(self.extx*self.exty*k+self.extx*j+i)
                            if (x0>=0 and x1<self.extx and y0>=0 and y1<self.exty and z0>=0 and z1<self.extz):
                                self.protein_phi[reorient(z0,y0,x0)] += weight* a * b * c
                                self.protein_phi[reorient(z1,y0,x0)] += weight * a * b * (1-c) 
                                self.protein_phi[reorient(z0,y1,x0)] += weight * a * (1-b) * c 
                                self.protein_phi[reorient(z0,y0,x1)] += weight * (1-a) * b * c 
                                self.protein_phi[reorient(z1,y1,x0)] += weight * a * (1-b) * (1-c)
                                self.protein_phi[reorient(z0,y1,x1)] += weight * (1-a) * (1-b) * c
                                self.protein_phi[reorient(z1,y0,x1)] += weight* (1-a) * b * (1-c)
                                self.protein_phi[reorient(z1,y1,x1)] += weight* (1-a) * (1-b) * (1-c)
                            else:
                                if (x0>=0 and x0<self.extx):
                                    if (y0>=0 and y0<self.exty): 
                                        if (z0>=0 and z0<self.extz):
                                            self.protein_phi[reorient(z0,y0,x0)] += weight* a * b * c
                                        if (z1>=0 and z1<self.extz):
                                            self.protein_phi[reorient(z1,y0,x0)] += weight * a * b * (1-c) 
                                    if (y1>=0 and y1<self.exty):
                                        if (z0>=0 and z0<self.extz): 
                                            self.protein_phi[reorient(z0,y1,x0)] += weight * a * (1-b) * c 
                                        if (z1>=0 and z1<self.extz):
                                            self.protein_phi[reorient(z1,y1,x0)] += weight * a * (1-b) * (1-c) 
                                if (x1>=0 and x1<self.extx):
                                    if (y0>=0 and y0<self.exty):
                                        if (z0>=0 and z0<self.extz):
                                            self.protein_phi[reorient(z0,y0,x1)] += weight * (1-a) * b * c 
                                        if (z1>=0 and z1<self.extz):
                                            self.protein_phi[reorient(z1,y0,x1)] += weight* (1-a) * b * (1-c)
                                    if (y1>=0 and y1<self.exty) :
                                        if (z0>=0 and z0<self.extz):
                                            self.protein_phi[reorient(z0,y1,x1)] += weight * (1-a) * (1-b) * c
                                        if (z1>=0 and z1<self.extz):
                                            self.protein_phi[reorient(z1,y1,x1)] += weight* (1-a) * (1-b) * (1-c)
                                            
    def __interpolate_structure(self,  structure,  volume):
        """
        This function interpolates translated rotated structure onto density map grid
        """
        weight = 1.0
        index = 0
        varp=0.0
        nvox = volume.extx * volume.exty * volume.extz
        structure.protein_phi = zeros(nvox)
        for coords in structure.to_file_coordinates:
            # compute position within grid - width will be taken from the CCP4 file
            gx = (coords[0]-volume.gridx)/volume.widthx
            gy = (coords[1]-volume.gridy)/volume.widthy
            gz = (coords[2]-volume.gridz)/volume.widthz
            x0 = floor(gx)
            y0 = floor(gy)
            z0 = floor(gz)
            x1 = x0+1
            y1 = y0+1
            z1 = z0+1
            #interpolate
            a = x1-gx
            b = y1-gy
            c = z1-gz
            weight = structure.atomic_weights[index]
            reorient = lambda k, j, i : int(volume.extx*volume.exty*k+volume.extx*j+i)
            if (x0>=0 and x1<volume.extx and y0>=0 and y1<volume.exty and z0>=0 and z1<volume.extz):
                structure.protein_phi[reorient(z0,y0,x0)] += weight* a * b * c
                varp+=weight * a * b * c * ((1-a)*(1-a)+(1-b)*(1-b)+(1-c)*(1-c))
                structure.protein_phi[reorient(z1,y0,x0)] += weight * a * b * (1-c) 
                varp+=weight * a * b * (1-c) * ((1-a)*(1-a)+(1-b)*(1-b)+c*c)
                structure.protein_phi[reorient(z0,y1,x0)] += weight * a * (1-b) * c 
                varp+=weight * a * (1-b) * c * ((1-a)*(1-a)+b*b+(1-c)*(1-c))
                structure.protein_phi[reorient(z0,y0,x1)] += weight * (1-a) * b * c 
                varp+=weight * (1-a) * b * c * (a*a+(1-b)*(1-b)+(1-c)*(1-c))
                structure.protein_phi[reorient(z1,y1,x0)] += weight * a * (1-b) * (1-c)
                varp+=weight * a * (1-b) * (1-c) * ((1-a)*(1-a)+b*b+c*c)
                structure.protein_phi[reorient(z0,y1,x1)] += weight * (1-a) * (1-b) * c
                varp+=weight * (1-a) * (1-b) * c * (a*a+b*b+(1-c)*(1-c))
                structure.protein_phi[reorient(z1,y0,x1)] += weight* (1-a) * b * (1-c)
                varp+=weight * (1-a) * b * (1-c) * (a*a+(1-b)*(1-b)+c*c)
                structure.protein_phi[reorient(z1,y1,x1)] += weight* (1-a) * (1-b) * (1-c)
                varp+=weight * (1-a) * (1-b) * (1-c) * (a*a+b*b+c*c)
            else:
                if (x0>=0 and x0<volume.extx):
                    if (y0>=0 and y0<volume.exty): 
                        if (z0>=0 and z0<volume.extz):
                            structure.protein_phi[reorient(z0,y0,x0)] += weight* a * b * c
                            varp+=weight * a * b * c * ((1-a)*(1-a)+(1-b)*(1-b)+(1-c)*(1-c))
                        if (z1>=0 and z1<volume.extz):
                            structure.protein_phi[reorient(z1,y0,x0)] += weight * a * b * (1-c)
                            varp+=weight * a * b * (1-c) * ((1-a)*(1-a)+(1-b)*(1-b)+c*c)
                    if (y1>=0 and y1<volume.exty):
                        if (z0>=0 and z0<volume.extz): 
                            structure.protein_phi[reorient(z0,y1,x0)] += weight * a * (1-b) * c
                            varp+=weight * a * (1-b) * c * ((1-a)*(1-a)+b*b+(1-c)*(1-c))
                        if (z1>=0 and z1<volume.extz):
                            structure.protein_phi[reorient(z1,y1,x0)] += weight * a * (1-b) * (1-c)
                            varp+=weight * a * (1-b) * (1-c) * ((1-a)*(1-a)+b*b+c*c)
                if (x1>=0 and x1<volume.extx):
                    if (y0>=0 and y0<volume.exty):
                        if (z0>=0 and z0<volume.extz):
                            structure.protein_phi[reorient(z0,y0,x1)] += weight * (1-a) * b * c 
                            varp+=weight * (1-a) * b * c * (a*a+(1-b)*(1-b)+(1-c)*(1-c))
                        if (z1>=0 and z1<volume.extz):
                            structure.protein_phi[reorient(z1,y0,x1)] += weight* (1-a) * b * (1-c)
                            varp+=weight * (1-a) * b * (1-c) * (a*a+(1-b)*(1-b)+c*c)
                    if (y1>=0 and y1<volume.exty) :
                        if (z0>=0 and z0<volume.extz):
                            structure.protein_phi[reorient(z0,y1,x1)] += weight * (1-a) * (1-b) * c
                            varp+=weight * (1-a) * (1-b) * c * (a*a+b*b+(1-c)*(1-c))
                        if (z1>=0 and z1<volume.extz):
                            structure.protein_phi[reorient(z1,y1,x1)] += weight* (1-a) * (1-b) * (1-c)
                            varp+=weight * (1-a) * (1-b) * (1-c) * (a*a+b*b+c*c)
            index+=1
        return varp
    
    def __interpolate_structure_blur(self,  structure,  volume):
        """
        This function interpolates translated rotated structure onto density map grid
        """
        weight = 1.0
        margin = 2
        minx, miny, minz, maxx, maxy, maxz=structure.return_structure_span()
        nvox = volume.extx * volume.exty * volume.extz
        structure.protein_phi = zeros(nvox)
        index = 0
        varp=0.0
        for coords in structure.to_file_coordinates:
            # compute position within grid - width will be taken from the CCP4 file
            gx = margin+(coords[0]-minx)/volume.width
            gy = margin+(coords[1]-miny)/volume.width
            gz = margin+(coords[2]-minz)/volume.width
            x0 = int(floor(gx))
            y0 = int(floor(gy))
            z0 = int(floor(gz))
            x1 = x0+1
            y1 = y0+1
            z1 = z0+1
            #interpolate
            a = x1-gx
            b = y1-gy
            c = z1-gz                            
            reorient = lambda k, j, i : int(volume.extx*volume.exty*k+volume.extx*j+i)
            weight = structure.atomic_weights[index]
            structure.protein_phi[reorient(z0,y0,x0)] += weight* a * b * c
            varp+=weight * a * b * c * ((1-a)*(1-a)+(1-b)*(1-b)+(1-c)*(1-c))
            structure.protein_phi[reorient(z1,y0,x0)] += weight * a * b * (1-c) 
            varp+=weight * a * b * (1-c) * ((1-a)*(1-a)+(1-b)*(1-b)+c*c)
            structure.protein_phi[reorient(z0,y1,x0)] += weight * a * (1-b) * c 
            varp+=weight * a * (1-b) * c * ((1-a)*(1-a)+b*b+(1-c)*(1-c))
            structure.protein_phi[reorient(z0,y0,x1)] += weight * (1-a) * b * c 
            varp+=weight * (1-a) * b * c * (a*a+(1-b)*(1-b)+(1-c)*(1-c))
            structure.protein_phi[reorient(z1,y1,x0)] += weight * a * (1-b) * (1-c)
            varp+=weight * a * (1-b) * (1-c) * ((1-a)*(1-a)+b*b+c*c)
            structure.protein_phi[reorient(z0,y1,x1)] += weight * (1-a) * (1-b) * c
            varp+=weight * (1-a) * (1-b) * c * (a*a+b*b+(1-c)*(1-c))
            structure.protein_phi[reorient(z1,y0,x1)] += weight* (1-a) * b * (1-c)
            varp+=weight * (1-a) * b * (1-c) * (a*a+(1-b)*(1-b)+c*c)
            structure.protein_phi[reorient(z1,y1,x1)] += weight* (1-a) * (1-b) * (1-c)
            varp+=weight * (1-a) * (1-b) * (1-c) * (a*a+b*b+c*c)
            index+=1
        return varp

    def calculate_cc_fast(self):
        #Generating protein map for the translated protein
        varp = self.__interpolate_structure(self.structure,  self.volume)
        varp = varp/self.structure.molecule_mass
        ##self.volume.write_to_volume_file(self.structure.protein_phi,  index)
        #If you want to have  smothed map
        if self.volume.resolution:
            self.volume.compute_Gaussian_kernel_acc(self.volume.resolution,  varp)
        return self.volume.calculate_correlation_fast(self.structure.protein_phi, self.structure.estimated_threshold)
    
    def blur_structure(self,  pdb_name):
        varp = self.__interpolate_structure_blur(self.structure,  self.volume)
        varp = varp/self.structure.molecule_mass
        #If you want to have  smothed map
        self.volume.compute_Gaussian_kernel_acc(self.volume.resolution, varp)
        simulated_protein_phi = self.volume.convolute_map_with_kernel_blur(self.structure.protein_phi)
        self.volume.write_to_volume_file(simulated_protein_phi,  pdb_name)
        return simulated_protein_phi
        
    def calculate_cc(self):
        """
        Taken almost entirely from UCSF chimera
        Arrays must be the same size!!!
        Numexpr might be used here
        """
        #Generating protein map for the translated protein
        self.__interpolate_structure(self.structure,  self.volume)
        ##self.volume.write_to_volume_file(self.structure.protein_phi,  index)
        #If you want to have  smothed map
        if self.volume.resolution:
            self.volume.compute_Gaussian_kernel(self.volume.resolution)
            ##self.structure.protein_phi = self.volume.convolute_map_with_kernel_fast(self.structure.protein_phi)
            self.structure.protein_phi = self.volume.convolute_map_with_kernel(self.structure.protein_phi)
            #Testing approach by Topf
            #self.structure.protein_phi=self.volume.compute_fourier_gaussian(self.structure.protein_phi,  self.volume.resolution,  0.225)
            #print self.structure.protein_phi
        ##self.volume.write_to_volume_file(self.structure.protein_phi,  1)
        threshold_indexes = self.structure.protein_phi>15.0
        protein_phi_above_thr = self.structure.protein_phi[threshold_indexes]
        map_phi_above_thr = self.volume.map_phi[threshold_indexes]
        overlap = dot(protein_phi_above_thr, map_phi_above_thr)
        norm1 = dot(protein_phi_above_thr, protein_phi_above_thr)
        norm2 = dot(map_phi_above_thr, map_phi_above_thr)
        n = len(protein_phi_above_thr)
        avg1 = sum(protein_phi_above_thr)/n
        avg2 = sum(map_phi_above_thr)/n
        d2 = (norm1 - n*avg1*avg1)*(norm2 - n*avg2*avg2)
        #Next few lines have been copied from UCSF chimera
        if d2 < 0:
            d2 = 0          # This only happens due to rounding error.
        d = sqrt(d2)
        if d == 0:
            cc = 0
        else:
            cc = (overlap - n*avg1*avg2)/d
        #print "DEBUG>>Chimera CC",  cc
        #situs_norm1 = norm(self.structure.protein_phi-avg1)
        #situs_norm2 = norm(self.volume.map_phi-avg2)
        #overlap = dot(self.structure.protein_phi-avg1, self.volume.map_phi-avg2)
        #cc = overlap/(situs_norm1*situs_norm2)
        return cc


    def __write_to_volume_file(self,  origx,  origy,  origz,  g_extx,  g_exty,  extz):
        """
        There is no actuall need for the dumping density into file but it might be hellpful for testing the algorithm.\
        Situs format at the moment. Testing purpose only
        """
        print "Writting density to file now"
        volume_file = open("volume_protein.sit",  "w")
        output_line ="%f %f %f %f %d %d %d\n" % (float(self.widthx),  float(origx), float(origy), float(origz), int(g_extx),  int(g_exty),  int(extz))
        volume_file.write(output_line)
        volume_file.write("\n")
        count = 0 
        for p in self.protein_phi:
            if (count+1)%10 == 0:
                output_line = " %10.6f \n" % float(p)
                volume_file.write(output_line)
            else:
                output_line = " %10.6f " % float(p)
                volume_file.write(output_line)
            count +=1
        volume_file.close()


if __name__ == "__main__":
    volume_file = sys.argv[2]
    pdb_file_name = sys.argv[1]
    threshold = float(sys.argv[3])
    resolution = float(sys.argv[4])
    volume = EMmap(volume_file,  threshold,  resolution)
    ##volume.extract_pseudo_atoms()
    volume.read_header_values()
    volume.read_volume_fast()
    volume.set_density_over_threshold()
    atomic_structure = AtomicStructure()
    atomic_structure.read(pdb_file_name)
    atomic_structure.estimated_threshold = threshold
    #subunit_pdb_file = "data/2G4C_AC.pdb"
    #This is the map threshold...
    corcoe = CorCoe(volume, atomic_structure)
    print corcoe.calculate_cc_fast()
