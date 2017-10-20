#!/usr/bin/python
"""
Library fo the reading procedures of the electron density maps formats
It is  partialy based on UCSF chimera and Situs code (some of the functions were imported) and there is no waranty for the full functionality
"""
__author__ = "Wojtek Potrzebowsk"
__copyright__ = "Copyright 2009"
__license__ = "GPL"
__version__ = "0.1.0"
__maintainer__ = "Wojtek Potrzebowski"
__email__ = "wpotrzeb@genesilico.edu.pl"
__status__ = "Prototype"

import os.path
from numpy import uint8, int8, int16, int32,  float32, dtype,  fromstring,  array

from numpy import reshape
from numpy import where
from numpy import zeros

class CCP4:

    def __init__(self, file_path):
        self.volume = open(file_path, 'rb')
        self.volume.seek(0,2)                              # go to end of file
        file_size = self.volume.tell()
        self.volume.seek(0,0)                              # go to beginning of file
        # Infer file byte order from column axis size nc.  Requires nc < 2**16
        # Was using mode value but 0 is allowed and does not determine byte order.
        self.swap_bytes = 0
        nc = self.read_values(self.volume, int32, 1)
        self.swap_bytes = not (nc > 0 and nc < 65536)
        self.volume.seek(0,0)
        self.v = self.read_header_values(self.volume, file_size)
        densities = []
        self.data_offset = self.volume.tell()
        esize = array((), float32).itemsize
        while (self.data_offset<file_size):
            ##densities.append(self.read_array(float32, 4,  self.swap_bytes))
            #This read might be buggy - have to keep an eye on it
            densities.extend(self.read_array(esize, 4,  self.swap_bytes))
            self.data_offset = self.volume.tell()
        self.volume.close()
        self.densities = array(densities)
		
        self.simulated_density_coords = {}
        self.density_coords = {}
        self.dotvalues = []
        self.width = 0.
        self.threshold = 0.0
        
    def read_header_values(self, file, file_size):
        MRC_USER = 29
        CCP4_USER = 15
        MRC_NUM_LABELS = 10
        MRC_LABEL_SIZE = 80
        MRC_HEADER_LENGTH = 1024
        i32 = int32
        f32 = float32
        v = {}
        v['nx'], v['ny'], v['nz'] = self.read_values(file, i32, 3)
        v['mode'] = self.read_values(file, i32, 1)
        v['mxstart'], v['mystart'], v['mzstart'] = self.read_values(file, i32, 3)
        v['mx'], v['my'], v['mz'] = self.read_values(file, i32, 3)
        v['xlen'], v['ylen'], v['zlen'] = self.read_values(file, f32, 3)
        v['alpha'], v['beta'], v['gamma'] = self.read_values(file, f32, 3)
        v['mapc'], v['mapr'], v['maps'] = self.read_values(file, i32, 3)
        v['amin'], v['amax'], v['amean'] = self.read_values(file, f32, 3)
        v['ispg'], v['nsymbt'] = self.read_values(file, i32, 2)
        v['lskflg'] = self.read_values(file, i32, 1)
        v['skwmat'] = self.read_values(file, f32, 9)
        v['skwtrn'] = self.read_values(file, f32, 3)
        v['user'] = self.read_values(file, i32, CCP4_USER)
        v['map'] = file.read(4)   # Should be 'MAP '.
        v['machst'] = self.read_values(file, i32, 1)
        v['rms'] = self.read_values(file, f32, 1)
        v['type'] = 'ccp4'
        v['nlabl'] = self.read_values(file, i32, 1)
        labels = []
        for i in range(MRC_NUM_LABELS):
            labels.append(file.read(MRC_LABEL_SIZE))
        v['labels'] = labels
        # Catch incorrect nsymbt value.
        if v['nsymbt'] < 0 or v['nsymbt'] + MRC_HEADER_LENGTH > file_size:
            raise SyntaxError, ('MRC header value nsymbt (%d) is invalid'
                          % v['nsymbt'])
        v['symop'] = file.read(v['nsymbt'])
        return v
    
    def read_values(self, file, etype, count):
        esize = array((), etype).itemsize
        string = file.read(esize * count)
        if len(string) < esize * count:
            raise SyntaxError, ('MRC file is truncated.  Failed reading %d values, type %s' % (count, etype.__name__))
        values = self.read_values_from_string(string, etype, count)
        return values

    def read_values_from_string(self, string, etype, count):
        values = fromstring(string, etype)
        if self.swap_bytes:
            values = values.byteswap()
        if count == 1:
            return values[0]
        return values

    def read_array(self, esize,  count,  swap):
        #It looks esize always equals 4
        string = self.volume.read(esize * count)
        currfloat = self.read_values_from_string(string, float32, count)
        #print currfloat
        #Swaping
        if swap == 1: 
            cptr = currfloat
            tmp = cptr[0]
            cptr[0]=cptr[3]
            cptr[3]=tmp
            tmp = cptr[1]
            cptr[1]=cptr[2]
            cptr[2]=tmp
            currfloat = cptr
        return currfloat
    
    def return_density_values(self):
        return self.densities
    
    def return_header_values(self):
        return self.v
				
    def read_simulated_volume_fast(self,  kmax, first_complex):
        """
        Reading volume and returning dictionary [x,y,z]=density
        """
        #TODO: In the mesher this values are reread. It can be done once
        v = self.return_header_values()
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
        extxy = self.extx * self.exty
		#return max(self.widthx, self.widthy, self.widthz)
		
        self.densities = reshape(self.densities, nvox)
        #############################
        molecule_mass = kmax*(first_complex.volume*0.81)
        Nv= float(molecule_mass)/(0.81*self.widthx*self.widthy*self.widthz)
        Nv = int(Nv*1.2)
        dot = self.densities[self.densities>0.0]
        dot = sorted(dot)
        self.threshold = round(dot[-Nv], 2)
        #structure.estimated_threshold = threshold
            
        print "complex MASS:", molecule_mass
        print "compelx VOLUME", first_complex.volume
        print "density threshold for %s volume is: %s"%(kmax, self.threshold) 
        
        #############################
        print "<<DEBUG: Threshold of simulated map is set to",  self.threshold
        counts = where(self.densities>self.threshold)
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
            density_value = self.densities[count]
            self.simulated_density_coords[(float(x), float(y), float(z))] = density_value
            self.dotvalues.append(density_value)
            count = count +1
            #print x,y,z, density_value
            
        return self.simulated_density_coords
			
    def get_map_threshold(self):
        return self.threshold
    
    def read_volume_fast(self, threshold):
        """
        Reading volume and returning dictionary [x,y,z]=density
        """
        v = self.return_header_values()
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
		
        densities = reshape(self.densities, nvox)
        counts = where(self.densities>threshold)
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
            density_value = self.densities[count]
            self.density_coords[(float(x), float(y), float(z))] = density_value
            #Needed for neighboor search I guess there should be interaface between EMmap and AtomicStrcucture  class
            #self.neighbor_list.append(Neighbor([float(x), float(y), float(z)]))
            self.dotvalues.append(density_value)
            #self.map_phi[count]+=float(density_value)
            count = count +1
            #print x,y,z, density_value
        #self.density_coords_bkp = deepcopy(self.density_coords)
        return self.density_coords
    
if __name__=="__main__":
    volume_name = "2REC.2.1.3.-15.1.1.1.mrc"
    ccp4 = CCP4(volume_name)
    densities = ccp4.return_density_values()
    vals = ccp4.read_simulated_volume_fast(volume_name, 1.6)
    
