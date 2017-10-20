# Copyright '2003-2006 Marta Kaczor, Michal J. Gajda, Anastasia Bakulina
# See license.txt

import PModel  #getting access to stride_dict
from Restraint import Restraint
import sys

class AccessRestraint(Restraint):
    'checking residue accessibility (%)'
    
    type=["Accessibility"]
    
    def __init__(self, res, lowbound, highbound, weight, id=''):
       #print "$$$$$$$$$$$$$$$$$$$$$$44", res
       self.res = res
       self.lowbound = lowbound
       self.highbound = highbound
       if self.lowbound > self.highbound: (self.lowbound, self.highbound) = (self.highbound, self.lowbound)
       self.weight = float(weight)
       self.f_punish = self.F_PUNISH_LIN # default
       self.id = id
       
    def check(self, pmodel):

### should be not here:
       #somehow check if dssp has been called:
       try:
          test = pmodel.model.get_list()[0].get_list()[0].dssp
       except AttributeError:
          pmodel.dssp()


       if self.res.chain_id == '': self.res.chain_id = pmodel.first_chain_id
       try:
          tmp_res = self.res.find1res(pmodel)
          acc = tmp_res.dssp
       except AttributeError:
          error_message = "*** Accessibility cannot be calculated. "
          if tmp_res != None:
             error_message += "Residue "+str(self.res)+" in "+pmodel.filename+" has no dssp value"
          else:
             error_message += "Residue "+str(self.res)+" not found in "+pmodel.filename   
          print >> sys.stderr, error_message  
          return 0
       if acc<self.highbound and acc>self.lowbound:
          return 0
       else:
          return self.calc_punish(min(abs(acc-self.highbound),abs(acc-self.lowbound))) * self.weight 
            


    def autoname (self):
       return  str(self.res)+' '+str(self.lowbound)+'-'+str(self.highbound)+'%'
       

    def pymol_color(self):
       """ Returns name of the color for this restraint. """    
       return "blue"
      
    def pymol_script(self):
       if self.lowbound==0:
         bounds_label="SA<=%d%%" % self.highbound
       elif self.highbound==100:
         bounds_label="SA=%d%%" % self.lowbound
       else:
         bounds_label="%d<=SA<=%d" % (self.lobound, self.highbound)
       return """
flag ignore, not (%(res)s), set
delete indicate
select %(name)s, (%(res)s)
color %(colorName)s, %(name)s
label (%(central_atom)s), "%(bounds_label)s"
set transparency=0.5, %(name)s
set two_sided_lighting=1
set backface_cull=0
show surface, %(name)s
flag ignore, not (%(res)s), clear
"""        % { "name"          : self.pymol_name(),
               "res"           : self.res.pymol_selector(),
               "lowbound"      : self.lowbound,
               "highbound"     : self.highbound,
               "bounds_label"  : bounds_label,
               "colorName"     : self.pymol_color(),
               "central_atom"  : self.res.central_atom().pymol_selector()
             }
       
       
