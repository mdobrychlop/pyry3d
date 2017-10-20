# Copyright 2003-2006 Marta Kaczor, Michal J. Gajda, Anastasia Bakulina
# see license.txt
import PModel  
from Restraint import Restraint
import sys

class SecStructureRestraint(Restraint):
    'checking correspondence of secondary structure'

    type=["Secondary Structure","Simple"]
    struct_type = {'helix':['H','G','I'],'strand':['E', 'B', 'b'],'coil':['C', 'T'],'flexible':['C', 'T','D']}
    def __init__(self, sstype, res_range, weight, id = ''):
       self.sstype = sstype # list of ss types
       self.res_range = res_range
       self.weight = weight 
       self.id = id
       
          

    def check(self, pmodel):
       wr = 0
     
       try: 
          test = pmodel.model.get_list()[0].get_list()[0].secondary_structure
       except AttributeError:
          pmodel.stride()
       #
       
       if self.res_range.chain_id == '': self.res_range.chain_id = pmodel.first_chain_id
     
       if self.sstype in self.struct_type.keys():
          ok_struct = self.struct_type[self.sstype]
       else:
          pass
          #print about error

       for res in self.res_range.find_res(pmodel):
          try: 
             sec_struct = res.secondary_structure
          except AttributeError:
#             print >> sys.stderr, "Cannot find secondary structure for residue ", res
             sec_struct = None   
          #we suppose that if residue does not exist and sstype is 'D' it is ok:
          if (sec_struct not in ok_struct) and not (res.get_resname()=='nonexistent' and 'D' in ok_struct): 
             print >> sys.stderr, "Secondary structure", sec_struct, "for residue", res, "is not in", ", ".join(ok_struct)
             wr += 1
       
       p = self.calc_punish((float(wr)/float(len(self.res_range))) * self.PENALTY_MAX) * float(self.weight)
       return p
       


    def autoname (self):
       return self.sstype+" "+str(self.res_range)
       
       
    def pymol_color(self):
       """ Returns name of the color for this restraint. """    
       sstype=self.get_SS_type()
       if sstype == 'H':
         return "yellow" # Helix
       elif sstype == 'E':
         return "red"    # Beta
       elif sstype == 'C':
         return "gray"   # loop
    
    def get_SS_type(self):
       sstype = self.sstype.lower()
       print >> sys.stderr, self.sstype
       if self.sstype.lower() == 'helix':
         return "H"   # Helix
       elif self.sstype.lower() == 'strand':
         return "E"   # Beta
       elif self.sstype.lower() in ['coil', 'loop', 'flexible']:
         return "C"   # loop
       else:
         print >>sys.stderr, "Unknown ss type %s" % repr(self.sstype)
    
    def pymol_script(self):
       sstype = self.get_SS_type()
       if sstype=='E':
         sstype='S'
       return """
         select %(name)s, (%(residues)s)
         alter /%(residues)s/, ss='%(sstype)s'
         color %(colorName)s, %(name)s
         show cartoon, %(name)s
       """ % { "name"       : self.pymol_name(),
               "residues"   : self.res_range.pymol_selector(),
               "sstype"     : sstype,
               "colorName"  : self.pymol_color()
             }
       
