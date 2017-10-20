# Copyright 2003-2006 Marta Kaczor, Michal J. Gajda, Anastasia Bakulina
# see license.txt
import PModel  #getting access to stride_dict
import re

global global_index
global_index=0 # Unique index of a restraint

def number_to_alpha(number):
    base = 26
    if number == 0:
        return "A"
    result = []
    while number>0:
        result.append(chr(ord('A')+(number % base)))
        number=int(number/base)
    result.reverse()
    return "".join(result)

class Restraint:
    '''Base class for all restraints.
       Instances of this class can be created for debugging'''
    type=[]
    name="fill short description of restraint type in a subclass"
    """o type -- List of strings categorizing restraint.
                 [category,subcategory] is used for summary columns and their labels
       o PENALTY_MAX -- Maximum penalty for a restraint.
       o PENALTY_NODATA -- Penalty for missing data.
       o F_PUNISH_LIN -- Linear increase in penalty with growing distance.
       o F_PUNISH_SQR -- Quadratic increase in penalty with growing distance. """
    PENALTY_MAX=20 # Distance deviation by more 20 angstroms can be counted as far too far...
    PENALTY_NODATA=PENALTY_MAX
    F_PUNISH_LIN = 1
    F_PUNISH_SQR = 2
    

    def set_punish(self,fp=F_PUNISH_LIN):
       self.f_punish = fp
      
    def calc_punish(self,p):
       p = float(p)
       if self.f_punish == self.F_PUNISH_LIN: return p
       elif self.f_punish == self.F_PUNISH_SQR: return p**2
       assert(False) # Should never happen!
           
    
    def check(self, pmodel):
        """ Return restraint penalty for a given model.
        
        o pmodel -- PModel object holding a model.
        
        o result -- returns a value between 0 and PENALTY_MAX,
                    0 means satisfied restraint
                    and 100 means that model is not even close
                    or doesn't contain necessary data.
        """
        return self.PENALTY_MAX
    
    def check_residues(self, pmodel):
        """ Check existence of residues mentioned in the restraint. """
        return '' #still not implemented
        
    def get_id(self):
       id = self.id[0:-1].strip() #remove ':' and spaces
       return id
        
    def __init__(self,id = ''):
       """ Create restraint with all it's parameters """
       self.id = id
       
    def name(self):
       """ Label for a column with results of the restraint """
       if self.id=="" or self.id is None:
         return self.__class__.__name__+" "+str(self.autoname())
       else:
         return self.id
     
    def __repr__(self):
       vars = filter(lambda name: (type(getattr(self, name))==str or
                                   type(getattr(self, name))==int or
                                   type(getattr(self, name))==float) and
                                  name[0]!='_',
                     dir(self))
       vars_descriptions = map(lambda name: "%s=%s"
                                          % (name, repr(getattr(self, name))),
                               vars)
       return "%s(%s)" % (self.__class__.__name__, ', '.join(vars_descriptions))
    
    def pymol_name(self):
       """ Returns a name that contains just alphabetic and '_' characters """
       string=re.sub("[^A-Za-z_]+", "_", self.name()[1])
       global global_index
       string=string.split("_")[0]+"_"+number_to_alpha(global_index)
       global_index=global_index+1 # WARNING: Nasty side-effect here!
       return string
    
    def pymol_color(self):
       """ Returns name of the color for this restraint. """
       return "blue"
    
    def pymol_script(self):
       """ Returns string that visualizes restraint in PyMol

          Returns script fragment that may be used
          as PyMol code to visualize a restraint.
       """
       return ""

