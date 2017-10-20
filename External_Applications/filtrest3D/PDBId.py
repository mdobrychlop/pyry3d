#!/usr/bin/python
# Copyright 2003-2006 Marta Kaczor, Michal J. Gajda, Anastasia Bakulina
# see license.txt

import sys, string
#from Bio.SCOP.Raf import to_one_letter_code
from Modules.Constans.Constans  import to_one_letter_code
from AATranslator import AATranslator

class PDBId:
   'class for searching atoms and residues in PModel'
   
   def __init__(self,res_name = "",res_num = None,chain_id = "",atom_name = ""):
      self.res_num = int(res_num)
      self.res_name = string.strip(res_name)
      self.chain_id = chain_id
      self.atom_name = string.strip(atom_name)
#      print 'res_num = ', self.res_num, 'res_name = ',self.res_name, 'chainID = ', self.chain_id, self.atom_name   

   def residues(self):
      return PDBId(self.res_name, self.res_num, self.chain_id, "")

   def central_atom(self):
      if self.atom_name!="" and self.atom_name!=None:
          return self
      else:
          return PDBId(self.res_name, self.res_num, self.chain_id, "CA")

   def __residue_wrong_name__(self,pmodel,ret_residue):
      warning_message = "Warning! Wrong residue name (residue '"+str(self.res_num)+"'"
      if self.chain_id != "": warning_message += '"'+self.chain_id+'"'
      warning_message += ' in ' + pmodel.filename
      warning_message += ' not "' + self.res_name + '"'
      if ret_residue != None: warning_message +=' but '+ret_residue.get_resname()+")."
      else: warning_message += ")."
      print >> sys.stderr, warning_message             
   
   
   def __residue_not_found__(self,pmodel):
     print >> sys.stderr,"Error! residue %s not found in %s!" % (str(self), pmodel.filename)   
     
   def __atom_not_found__(self,pmodel):
     print >> sys.stderr,"Error! atom %s not found in %s!" % (str(self), pmodel.filename)

   def find1res(self,pmodel):
      'returns one residue - for assigning sec.structure, accessibility and so on'
      ret_residue = None
      for chain in pmodel.get_iterator():
         if chain.get_id() == self.chain_id or self.chain_id == '':
            for residue in chain.get_iterator():
               if (residue.get_id()[1] == self.res_num) and \
               ((string.strip(residue.get_resname()) == self.res_name
                or self.res_name) == ""):
                  return residue
      
      if len(self.res_name) == 1: # try to convert residue.get_resname() to one-letter code
         for chain in pmodel.get_iterator():
            if chain.get_id() == self.chain_id or self.chain_id == '':
               for residue in chain.get_iterator():
                  resname = string.strip(residue.get_resname())
                  if to_one_letter_code.has_key(resname):
                     residue_resname = to_one_letter_code[resname]
                     if (residue.get_id()[1] == self.res_num) and \
                            (string.strip(residue_resname) == string.strip(self.res_name)):
                        return residue

      #maybe wrong residue name, so search again without name check:
      for chain in pmodel.get_iterator():
         if chain.get_id() == self.chain_id or self.chain_id == '':
            for residue in chain.get_iterator():
               #print residue.get_id()[1], self.res_num
               if residue.get_id()[1] == self.res_num:
                  ret_residue = residue
                  self.__residue_wrong_name__(pmodel,ret_residue)
                  return ret_residue
                  
      #if nothing has been found:
      self.__residue_not_found__(pmodel)
      return None
   
            
      
   def find_res(self,pmodel):
      'returns list of residues' 
      chain_iterator = pmodel.get_iterator()
      list = []
      resnumlist = []
      for chain in pmodel.get_iterator():
         if chain.get_id() == self.chain_id or self.chain_id == '':
            for residue in chain.get_iterator():
               if ((residue.get_id()[1] == self.res_num)  or  (self.res_num == None)) and \
                ((string.strip(residue.get_resname()) == self.res_name) or (self.res_name == "")):
                  list.append(residue)    
      if len(list) > 0:           
         return list
      
      elif len(self.res_name) == 1: # try to convert residue.get_resname() to one-letter code
         for chain in pmodel.get_iterator():
            if chain.get_id() == self.chain_id or self.chain_id == '':
               for residue in chain.get_iterator():
                  if to_one_letter_code.has_key(string.strip(residue.get_resname())):
                     residue_resname = to_one_letter_code[string.strip(residue.get_resname())]
                     if ((residue.get_id()[1] == self.res_num)  or 
                         (self.res_num == None)) and \
                         ((string.strip(residue_resname) == self.res_name) or
                         (self.res_name == "")):
                        list.append(residue)    
         if len(list) > 0:
            return list
            
      #maybe wrong residue name, so search again without name check:
      for chain in pmodel.get_iterator():
         if chain.get_id() == self.chain_id or self.chain_id == '':
            for residue in chain.get_iterator():
               if ((residue.get_id()[1] == self.res_num)  or  (self.res_num == None)):
                  list.append(residue)    
      if len(list) > 0:           
         if len(list) == 1: self.__residue_wrong_name__(pmodel,list[0])
         else: self.__residue_wrong_name__(pmodel,None)
         return list
      else:
         self.__residue_not_found__(pmodel)
         return None         
      
        
   def find_at(self,pmodel):
      'returns list of atoms' 
      chain_iterator = pmodel.get_iterator()
      list = []
      resnumlist = []
      for chain in pmodel.get_iterator():
         if chain.get_id() == self.chain_id or self.chain_id == '':
            for residue in chain.get_iterator():
               for atom in residue.get_iterator():
                  if ((residue.get_id()[1] == self.res_num) or (self.res_num == None)) and \
                     ((string.strip(residue.get_resname()) == self.res_name) or (self.res_name == "")) and \
                     ((string.strip(atom.get_id()) == self.atom_name) or (self.atom_name == "")):
                        list.append(atom)    
                 
      if len(list) > 0:           
         return list
      
      elif len(self.res_name) == 1: # try to convert residue.get_resname() to one-letter code
         for chain in pmodel.get_iterator():
            if chain.get_id() == self.chain_id or self.chain_id == '':
               for residue in chain.get_iterator():
                  if to_one_letter_code.has_key(residue.get_resname()):
                     residue_resname = to_one_letter_code[residue.get_resname()]
                     for atom in residue.get_iterator():
                        if ((residue.get_id()[1] == self.res_num) or (self.res_num == None)) and \
                           ((string.strip(residue_resname) == self.res_name) or (self.res_name == "")) and \
                           ((string.strip(atom.get_id()) == self.atom_name) or (self.atom_name == "")):
                              list.append(atom)    

         if len(list) > 0:
            return list
      
      #maybe wrong residue name, so search again without name check:
      for chain in pmodel.get_iterator():
         if chain.get_id() == self.chain_id or self.chain_id == '':
            for residue in chain.get_iterator():
               for atom in residue.get_iterator():
                  if ((residue.get_id()[1] == self.res_num) or (self.res_num == None)) and \
                     ((string.strip(atom.get_id()) == self.atom_name) or (self.atom_name == "")):
                        list.append(atom)    
                        ret_residue = residue
      if len(list) > 0:
         if len(list) == 1: self.__residue_wrong_name__(pmodel,ret_residue)
         else: self.__residue_wrong_name__(pmodel,None)
         return list
      else:
         self.__atom_not_found__(pmodel)
         return None      
   
   def pymol_selector(self):
      "Generate PyMol selector for residue(s)"   
      result=[]
      if self.res_num!=None:
        result+=["resi %d" % self.res_num]
      if self.atom_name!=None and self.atom_name!="":
        result+=["name %s" % self.atom_name]
      if self.res_name!=None and self.res_name!="":
        name=self.res_name
        if len(name)==1:
          name=AATranslator().get3letter(name)
        result+=["resn %s" % name]
      if self.chain_id!=None and self.chain_id!="":
        result+=["chain %s" % self.chain_id]
      return " and ".join(result)
                      
    
   def __str__(self):
      
      if self.res_num == None: res_num = ''
      else: res_num = str(self.res_num)
      ret_str = self.res_name+res_num
      if (self.chain_id):
         ret_str += '\"'+self.chain_id+'\"'
      if (self.atom_name):
         ret_str += ' '+self.atom_name
      return ret_str 
    
   def __repr__(self):
      return "PDBId(%s)" % ', '.join(map(repr, [self.res_name, self.res_num, self.chain_id, self.atom_name]))

 
