# Copyright 2003-2006 Marta Kaczor, Michal J. Gajda
# see license.txt
import PModel  #getting access to stride_dict
from Restraint import Restraint
import sys

class OrRestraint(Restraint):
    'maxumum penalty of several restraints'
    
    type=["Complex","Or"]
    
    def __init__(self, restr_list,id = ''):
       self.restr_list = restr_list
       self.id = id
       
    def check(self, pmodel):

       check_results = []
       for restr in self.restr_list:
          check_results.append(restr.check(pmodel))
       return min(check_results)

    def set_punish(self, p):
      for r in self.restr_list:
        r.set_punish(p)

    def name (self):
       return ' or '.join(map(lambda r: str(r.name()),
                              self.restr_list))
       
    def autoname (self):
       return '_or_'.join(map(lambda r: str(r.autoname()),
                              self.restr_list))
