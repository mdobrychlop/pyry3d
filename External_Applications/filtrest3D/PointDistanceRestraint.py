# Copyright '2011 by Joanna M. Kasprzak
# See license.txt

import PModel  #getting access to stride_dict
from Restraint import Restraint
import sys
from PDBId import PDBId



class PointDistance(Restraint):
    'check distance between two sets (or pair) of residues'
    type=["PointDistnace"]
    
    def __init__(self, reslist1, point, relation, distance, weight,id=''):
       self.reslist1 = reslist1 #list of PDBId
       self.point    = point
       self.dist     = distance
       self.relation = relation
       self.weight   = weight 
       self.id       = id
       
    def __str__(self):
      return "%s %s %s %s %s %s" % \
              (self.reslist1, self.point, self.dist, self.relation, self.weight, self.id)
       
    def autoname (self):
       return self.__reslist2str__(self.reslist1)+'-'+str(self.point)+" "+self.relation+str(self.dist)+'w = '+str(self.weight)
         
    def __reslist2str__ (self,reslist):
       str1 = '('
       for res in reslist:
          str1 += str(res)+','
       str1 += ')'   
       return str1