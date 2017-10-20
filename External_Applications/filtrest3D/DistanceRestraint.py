# Copyright '2006 Marta Kaczor, Michal J. Gajda, Anastasia Bakulina
# See license.txt

import PModel  #getting access to stride_dict
from Restraint import Restraint
import sys
from PDBId import PDBId

def list_middle_item(l):
  index=len(l)/2
  return l[index]

class DistanceRestraint(Restraint):
    'check distance between two sets (or pair) of residues'
    type=["Distance"]
    
    def __init__(self, reslist1, reslist2, relation, distance, weight,id=''):
       self.reslist1 = reslist1 #list of PDBId
       self.reslist2 = reslist2
       self.dist = distance
       self.relation = relation
       self.weight = weight 
       self.id = id
         
    def __check_pair__(self,coord1,coord2,dist,comp):
# return punish value without f_punish   
       (x1,y1,z1) = coord1
       (x2,y2,z2) = coord2
       s = ((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2) ** 0.5
       dif = 0
       if comp == '<' or comp == '<=':      dif = s - dist
       elif comp == '>' or comp == '>=':    dif = dist - s
       elif comp == '=':                    dif = abs (dist - s)   
       else: print >> sys.stderr,'*** Something really wrong with the program! This point should never be reached'
       
       if (dif > 0): 
          return dif
       else:
          return 0   
        
    
    def __check_res__(self,pmodel,res,at):
       return 0   


    def check(self, pmodel):
       punish = []
       punish.append(self.PENALTY_MAX)       
   
       try:    
         try:
          for res_id1 in self.reslist1:
             for at1 in res_id1.find_at(pmodel):
                for res_id2 in self.reslist2:
                   for at2 in res_id2.find_at(pmodel):
                      p = self.__check_pair__(at1.get_coord(), at2.get_coord(), self.dist, self.relation)
                      p = self.calc_punish(p)
                      punish.append(p)
       
          score = min(punish)   
          return score    
         except TypeError:
          print >> sys.stderr, "*** Distance cannot be calculated."
          return 0
       except AttributeError: 
         print map(lambda l:l.__class__, [res_id1, at1, res_id2, at2, self.reslist1, self.reslist2])
         raise
#    def check_residues(self, pmodel):
#       pass
#       # return warning


    def pymol_color(self):
        return "violet"
    
    def pymol_script(self):
       # TODO: fix bad treatment of residue PDBId
       return """
         select %(name)s_residues = ((%(residues1)s) or (%(residues2)s))
         #show sticks, %(name)s_residues
         select %(name)s_atoms, ((%(target1)s) or (%(target2)s))
         show spheres, %(name)s_atoms
         distance %(name)s = (%(target1)s),(%(target2)s)
         color purple, %(name)s_atoms
         color red, %(name)s_residues
         color %(colorName)s, %(name)s 
       """ % { "name"        : self.pymol_name(),
               "residues1"   : "("+(') or ('.join(map(lambda p: p.residues().pymol_selector(),
                                                      self.reslist1)))+")",
               "residues2"   : "("+(') or ('.join(map(lambda p: p.residues().pymol_selector(),
                                                      self.reslist2)))+")",
               "target1"     : list_middle_item(self.reslist1).central_atom().pymol_selector(),
               "target2"     : list_middle_item(self.reslist2).central_atom().pymol_selector(),
               "colorName"   : self.pymol_color()
             }
             
       
       

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
       return self.__reslist2str__(self.reslist1)+'-'+self.__reslist2str__(self.reslist2)+self.relation+str(self.dist)+'w = '+str(self.weight)

