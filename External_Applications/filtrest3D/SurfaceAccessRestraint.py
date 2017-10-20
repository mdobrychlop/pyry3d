# Copyright '2006 Marta Kaczor, Michal J. Gajda, Anastasia Bakulina
# See license.txt

import PModel  #getting access to stride_dict
from Restraint import Restraint
import sys
from PDBId import PDBId



class SurfaceAccess(Restraint):
    'check distance between two sets (or pair) of residues'
    type=["SurfaceAccess"]
    
    def __init__(self, reslist1, relation, distance, weight,id=''):
       self.reslist1 = reslist1 #list of PDBId
       self.dist     = distance
       self.relation = relation
       self.weight   = weight 
       self.id       = id
       
    def __str__(self):
      return "%s %s %s %s %s" % \
              (self.reslist1, self.dist, self.relation, self.weight, self.id)
       
    def autoname (self):
       return self.__reslist2str__(self.reslist1)+'-'+self.relation+str(self.dist)+'w = '+str(self.weight)
         
    def __reslist2str__ (self,reslist):
       str1 = '('
       for res in reslist:
          str1 += str(res)+','
       str1 += ')'   
       return str1