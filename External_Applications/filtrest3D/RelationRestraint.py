# Copyright '2006 Marta Kaczor, Michal J. Gajda, Anastasia Bakulina
# See license.txt

import PModel  #getting access to stride_dict
from Restraint import Restraint
import sys
from PDBId import PDBId

def list_middle_item(l):
  index=len(l)/2
  return l[index]

class Relation(Restraint):
    'checks whether relation between distances between two sets (or pair) of residues is fulfilled'
    type=["Distance"]
    
    def __init__(self, reslist1, reslist2, reslist3, reslist4, relation, weight,id=''):
       self.reslist1 = reslist1 #list of PDBId
       self.reslist2 = reslist2
       self.reslist3 = reslist3
       self.reslist4 = reslist4
       self.relation = relation
       self.weight = weight 
       self.id = id
       
    def __str__(self):
      return "%s %s %s %s %s %s %s" % \
              (self.reslist1, self.restlist2, self.restlist3, self.restlist4, self.relation, self.weight, self.id)
            
    def __reslist2str__ (self,reslist):
       str1 = '('
       for res in reslist:
          str1 += str(res)+','
       str1 += ')'   
       return str1

    def short_name(self):
       if self.id!="":
         return self.id
       else:
         return self.name()[1:]

    def autoname (self):
       return self.__reslist2str__(self.reslist1)+'-'+self.__reslist2str__(self.reslist2) #+\
                                  #(self.reslist3)+'-'+self.__reslist2str__(self.reslist4)+str(self.relation)+'w = '+str(self.weight)

