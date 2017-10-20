#!/usr/bin/env python
# -*- coding: utf-8 -*-

#biopython
from Bio import PDB
from Bio.PDB import PDBParser
from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure
#standard library
import glob, sys, os
#PyRy3D modules
from Modules.Error.Errors       import ShapeError
from Pyry_cleanPDB              import run_cleanPDB

#class ShapeError(Exception): pass

class SAXSShape(object):
    
    def __init__(self):
        self.da_nr = 0                     #number of dummy atoms
        self.da_radius = 0.0               #radius of a single dummy atom
        self.gradius = 0.0                 #radius of gyration
        self.shape_file = None             #saxs shape file name (cleaned)
        self.dummy_atoms = []
        
    def clean_saxs_file(self,filename):
        """
        """
        #clean dammin/f file
        self.shape_file = str(filename)+".pyry"
        #os.system("python cleanPDB.py "+str(filename)+" "+self.shape_file)
        run_cleanPDB(str(filename), self.shape_file)
        
    def extract_data(self, saxs_file):
        """
        """
        fh = open(saxs_file, "r")
        
        for line in fh:
            if line.startswith("REMARK"):
                l = line.split()
                if ("output" in l) and ("phase" in l):
                    self.da_nr = int(l[8])
                elif ("radius" in l) and ("gyration" in l):
                    self.gradius = float(l[7])
            if line.startswith("ATOM"):
                l = line.split()
                at = DummyAtom()
                #at.coord = [float(l[5]), float(l[6]), float(l[7])]
                at.coord = [float(line[30:37]), float(line[38:45]), float(line[46:53])]
                at.number = line[6:10]
                self.dummy_atoms.append(at)
        #if int(self.da_nr) != len(list(self.dummy_atoms)):
        #    raise ShapeError("Different number of dummy atoms in saxsfile and dummy atoms list")

        #if self.da_radius <= 0.0: raise ShapeError("File with ab initio reconstruction does not have information about the radius of dummy atoms")
        fh.close()
        
    
    def get_data_from_dammif_file(self, filename):
        """
        """
        self.clean_saxs_file(filename)
        #clean dammin/f file
        #out_saxs_filename = str(filename)+".pyry"
        
        #run_cleanPDB(str(filename), self.shape_file)
        #os.system("python cleanPDB.py "+str(filename)+" "+self.shape_file)
        
        #self.dummy_atoms = PDBParser().get_structure("saxs_shape", self.shape_file)
        
        
        self.extract_data(self.shape_file)
        
class DummyAtom:
    
    def __init__(self):
        self.coord = []
        self.number = None

if __name__=='__main__':
    
    #get mapfile
    saxs_file = sys.argv[1]
    
    
    #retrieve data from ccp4 file
    mapcomponent = Shape()
    mapcomponent.get_data_from_saxs_file(saxs_file)
    #radius = mapcomponent.mapgrid.radius 
    #st = Statistic(mapcomponent.mappoints, 1)
    #thres = st.generate_volume_chart(mapcomponent, cmplxv, radius)
    #print "Threshold", thres
    
    
    
