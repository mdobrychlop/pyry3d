#!/usr/bin/python
# Copyright 2003-2006 Anastasia Bakulina, Michal J. Gajda
# see license.txt

from PDBId import PDBId
import sys, string
from Bio.PDB.Atom import Atom
#from genesilico.structure.SMCRA_extended import Atom

class PDBrange(PDBId):
   'class for searching atoms and residues in PModel'
   
   def __init__(self,res1_name = "",res1_num = None,res2_name = "",res2_num = None,chain_id = "",atom_name = ""):
      self.res1_num = int(res1_num)   
      self.res1_name = string.strip(res1_name)
      self.res2_num = int(res2_num)
      self.res2_name = string.strip(res2_name)
      if self.res2_num < self.res1_num: 
         (self.res1_num, self.res2_num) = (self.res2_num, self.res1_num)
      self.chain_id = chain_id
      self.atom_name = string.strip(atom_name)
      #residue names are ignored

#   def __resrange_not_found__(self,pmodel):
#      print >> sys.stderr, "Error! No residues from range "+str(self)+" was found in "+pmodel.filename+"."
   
   def residues(self):
      return PDBrange(self.res1_name, self.res1_num, self.res2_name, self.res2_num, self.chain_id, "")
  
   def central_atom(self):
      resnum=self.res1_num+len(self)/2
      return PDBId(res_num=resnum, chain_id=self.chain_id).central_atom()
  
   def find_res(self,pmodel):
      'return list of residues' 
      if self.res1_num == None or self.res2_num == None: return None 
      list = []
      for chain in pmodel.get_iterator():
         if (chain.get_id() == self.chain_id) or (self.chain_id == None) or (self.chain_id == ''): 
            for residue in chain.get_iterator():
               for res_num in range(self.res1_num, self.res2_num+1):
                  if residue.get_id()[1] == res_num:
                     list.append(residue)    
      if len(list) == 0:
         print >> sys.stderr, "Error! No residues from range "+str(self)+" was found in "+pmodel.filename+"."
      return list
   
   def find_at(self,pmodel):
      'return list of atoms' 
      list = []
      for chain in pmodel.get_iterator():
         if (chain.get_id() == self.chain_id) or (self.chain_id == None) or (self.chain_id == ""):
            for residue in chain.get_iterator():
               for res_num in range(self.res1_num, self.res2_num+1):
                  if (residue.get_id()[1] == res_num):
                      for atom in residue.get_iterator():
                         if (string.strip(atom.get_id()) == self.atom_name):  
                           
                            list.append(atom)
                            if atom.__class__!=Atom:
                               print atom, type(atom), atom.__class__, residue, type(residue)

      if len(list) == 0:
         print >> sys.stderr, "Error! No atoms from range "+str(self)+" was found in "+pmodel.filename+"."
      return list
   
   def __str__(self):
      if self.res1_num == None: res1_num = ''
      else: res1_num = str(self.res1_num)
      if self.res2_num == None: res2_num = '' 
      else: res2_num = str(self.res2_num)
      ret_str = self.res1_name+res1_num+'-'+self.res2_name+res2_num   
      if (self.chain_id):
         ret_str += '\"'+self.chain_id+'\"'
      if (self.atom_name):
         ret_str += ' '+self.atom_name
      return ret_str 
    
   def __repr__(self):
      #return "PDBId(%s)" % ', '.join(map(repr, [res1_name, res1_num, res2_name, res2_num, chain_id, atom_name]))
      return "PDBId(%s)" % ', '.join(map(repr, [self.res1_name,self.res1_num, self.res2_name, self.res2_num, self.chain_id, self.atom_name]))
   
   def pymol_selector(self):
      "Generate PyMol selector for residue(s)"   
      result=[]
      if self.res1_num!=None or self.res2_num!=None:
        resinums="resi "
        if self.res1_num!=None:
          resinums+=str(self.res1_num)
        resinums+="-"
        if self.res2_num!=None:
          resinums+=str(self.res2_num)
        result+=[resinums]
      if self.atom_name!=None and self.atom_name!="":
        result+=["name %s" % self.atom_name]
      if self.chain_id!=None and self.chain_id!="":
        result+=["chain %s" % self.chain_id]
      return " and ".join(result)
    
   def __len__(self):
      return self.res2_num - self.res1_num +1
