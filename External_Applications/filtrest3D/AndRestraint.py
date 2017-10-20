# Copyright '2006 Marta Kaczor, Michal J. Gajda, Anastasia Bakulina
# See license.txt

import PModel  #getting access to stride_dict
from Restraint import Restraint
import sys

class AndRestraint(Restraint):
    'maxumum penalty of several restraints'
    
    type=["Complex","And"]
    
    def __init__(self, restr_list, id = ''):
       self.restr_list = restr_list
       self.id = id
       
    def check(self, pmodel):

       check_results = []
       for restr in self.restr_list:
          check_results.append(restr.check(pmodel))
       return max(check_results)

    def set_punish(self, p):
      for r in self.restr_list:
        r.set_punish(p)

    def name (self):
       return " and ".join(map(lambda r: str(r.name()), self.restr_list))
       
    def name (self):
       return "_and_".join(map(lambda r: str(r.name()), self.restr_list))
       
       
