# Copyright 2003-2006 Marta Kaczor, Michal J. Gajda, Anastasia Bakulina
# see license.txt
# -----------------------------------------------------------------------------
# RestraintsParser.py
#
# 
# Parser do plikow z wiezami
# -----------------------------------------------------------------------------
### new version of parser, returns list of Restraint objects
### Restraint objects are created with PDBId objects

# raz uzh tut tak mnogo kommetariev na pol'skom, dobavlju odin russkij, pust' potom dumajut chto tut napisano :)

import lex
import yacc
import re
#from Bio.SCOP.Raf import to_one_letter_code
import sys, os

from PDBId                  import PDBId
from PDBrange               import PDBrange

from DistanceRestraint      import DistanceRestraint
from Restraint              import Restraint
from SecStructureRestraint  import SecStructureRestraint
from AccessRestraint        import AccessRestraint
from AndRestraint           import AndRestraint
from OrRestraint            import OrRestraint
from Edm3SOMRestraint       import Edm3SOMRestraint
from EdmSitusRestraint      import EdmSitusRestraint
from SSPercentageRestraint  import SSPercentageRestraint
from SurfaceAccessRestraint import SurfaceAccess
from PointDistanceRestraint import PointDistance
from SymmetryRestraint      import Symmetry

class FormatError(Exception): pass

class Parser:

   tokens = ()
   reserved = {}
   
   def __init__(self):
      lex.lex(module=self)
      yacc.yacc(module=self)
      
   def run(self, infile):
      f = open(infile, 'r') 
      text = f.read()
      f.close()
      yacc.parse(text)


class RestraintsParser(Parser):
   def run(self, infile):
      f = open(infile, 'r') 
      text = f.read()
      f.close()
      if len(text)==0:
         yacc.parse(' ')  # bo dla pustego pliku jest blad: "No input string given with input()"
      else: yacc.parse(text)
      return (self.names, self.restrlist, self.symrestrlist)		
	
   tokens = (
      'ID', 'INT', 'FLOAT', 'STRING', 'CHAIN_ID',
      'LPAREN', 'RPAREN', 'COMMA',
      'DASH','WEIGHT', 
      'DEFAULT_WEIGHT',
      'ASSIGN',
      'COMP_LT', 'COMP_LE', 'COMP_EQ', 'COMP_GE', 'COMP_GT',
      'HELIX', 'STRAND', 'COIL','FLEXIBLE', 'DIST', 'SA', 'PD', 'SYM'
      'ACCESS','AND', 'OR','RID','EDM', 'SS_COMPOSITION' 
      )
	
   reserved = {
      'helix' : 'HELIX',
      'strand' : 'STRAND',
      'coil' : 'COIL',
      'flexible' : 'FLEXIBLE',
      'dist' : 'DIST',
      'pointdist' : 'PD',
      'surface_access' : 'SA',
      'symmetry' : 'SYM',
      'weight' : 'WEIGHT',
      'DEFAULT_WEIGHT' : 'DEFAULT_WEIGHT',
      'access' : 'ACCESS',
      'and': 'AND',
      'or':'OR',
      'edm':'EDM',
      'ss_composition': 'SS_COMPOSITION'
      }

   	
   t_LPAREN = r'\('
   t_RPAREN = r'\)'
   t_COMMA  = r'\,'
   t_DASH   = r'-'
   t_ASSIGN = r'<-'
   
   
   t_COMP_LT = r'<'
   t_COMP_LE = r'<='
   t_COMP_EQ = r'='
   t_COMP_GE = r'>='
   t_COMP_GT = r'>'
   
#   t_SIGN_ANY = '\*'
   
   t_ignore = ' \t'


      

   def t_STRING(self, t):
      r"'[^'\n]+'"
      try: t.value = t.value[1:-1]
      except ValueError:
         print "Problems with string?", t.value
         t.value = ""
      return t
    
   def t_FLOAT(self, t):
      r'\d+\.\d+'
      try: t.value = float(t.value)
      except ValueError:
         print "Float value too large", t.value
         t.value = 0.0
      return t


   def t_INT(self, t): 
      r'\d+'
      
      try: t.value = int(t.value)
      except ValueError:
         print "Integer value too large", t.value
         t.value = 0
      return t


  
   def t_RID(self, t):
      r'[a-zA-Z][a-zA-Z0-9_*]*\s*:' #now it is OK
      t.type = self.reserved.get(t.value, 'RID')
      t.value = t.value[:-1]
      return t # chop colon




   def t_ID(self, t):
#      r'[a-zA-Z_]+'
      r'(?:\#[a-zA-Z0-9*_]+\#) | (?:\#[a-zA-Z0-9\'_]+\#)|(?:[a-zA-Z_]+)'
      if t.value[0]=='#':
        t.value=t.value[1:]
      if t.value[-1]=='#':
        t.value=t.value[:-1]
      t.type = self.reserved.get(t.value, 'ID')
      return t

   def t_CHAIN_ID(self, t):
      r'\"[a-zA-Z0-9_]\"'
      t.value=t.value[1]
      return t

##############################


   def t_newline(self, t):
      r'(?:\n|\r)+'
      t.lineno += t.value.count('\n')


   def t_comment_1_line(self, t):
      r'//.*(?:\n|\r)'
      t.lineno += t.value.count('\n')

   def t_comment_multi_line(self, t):
      r'\(\*(.|\n|\r)*?\*\)'
      t.lineno += t.value.count('\n')

   def t_error(self, t):
      print "Illegal character %s" % repr(t.value[0])
      t.skip(1)

    
   names = {}
   restrlist = []
   DEFAULT_WEIGHT = { 'w' : 1}



   def p_main_0(self, p):
      'main : '

   def p_main_def_w(self, p):
      'main : def_weight main'

   def p_def_weight(self, p):
      'def_weight : DEFAULT_WEIGHT ASSIGN number'
      self.DEFAULT_WEIGHT['w'] = p[3]

   def p_main_1(self, p):
      'main : restr main'
      self.restrlist += p[1]

   def p_restr(self,p):
      '''restr : secstruct 
               | dist
               | surface_access
               | pointdist
               | access 
               | and
               | or
               | edm
               | ss_composition'''
      p[0] = p[1]                   


######################### EDM restraints

   edm_classes = {'3SOM':Edm3SOMRestraint,'Situs':EdmSitusRestraint}
   
   def p_edm_restr(self,p):
      'edm : EDM LPAREN opt_restrid STRING STRING RPAREN'
      edm_class = self.edm_classes.get(p[4],None)
      if edm_class == None:
         print '*** Unknown EDM docking method:', p[4]
         p[0] = []
      else:   
         restr = edm_class(p[5],p[3])
#         self.restrlist.append(restr)
         p[0] = [restr]

######################### Logical restraints
   
   def p_list_restr(self,p):
      'list_restr : restr list_restr'
      p[0] = p[1]+p[2]
      
   def p_list_restr_empty(self,p):   
      'list_restr :'
      p[0] = []
      
   def p_and(self,p):
      'and : AND opt_restrid LPAREN list_restr RPAREN' 
      restr = AndRestraint(p[4],p[2])
#      self.restrlist.append(restr)
      p[0] = [restr]

   def p_or(self,p):
      'or : OR opt_restrid LPAREN list_restr RPAREN' 
      restr = OrRestraint(p[4],p[2])
#      self.restrlist.append(restr)
      p[0] = [restr]
      

############################ Secondary structure restraints

   def p_secstruct(self,p):
      'secstruct : secstruct_type LPAREN list_el RPAREN'
      p[0] = []
      for (id,res_range, chain, weight) in p[3]:
         restr = SecStructureRestraint(p[1], res_range, weight,id)
#         self.restrlist.append(restr)
         p[0].append(restr)
         
   def p_secstruct_type(self,p):
      '''secstruct_type : HELIX 
                        | STRAND
                        | COIL
                        | FLEXIBLE'''
      p[0] = p[1]                 
   

   def p_list_el_0(self, p):
      'list_el : '
      p[0] = []
      
   def p_list_el_1(self, p):
      'list_el : el list_el'
      p[0] = p[1]+p[2]   


   def p_el_ch(self, p):         
      'el : opt_restrid res DASH res opt_chain opt_weight'     
      res_range = PDBrange(p[2][0],p[2][1],p[4][0],p[4][1],p[5])
      p[0] = [(p[1],res_range,p[5],p[6])]

########################weight, can be used in any restraint

   def p_opt_weight_empty(self,p):
      'opt_weight : '
      p[0] = self.DEFAULT_WEIGHT['w']

   def p_opt_weight(self,p):
      'opt_weight : WEIGHT ASSIGN number'
      p[0] = p[3]    

#######################Distance and fuzzy distance restraints

   def p_dist(self, p):
      'dist : DIST LPAREN dist_list RPAREN'
      p[0] = p[3]

   def p_dist_list_0(self, p):
      'dist_list : '
      p[0] = []       
  
   def p_dist_list_1(self, p):
      'dist_list : dist_el dist_list'
      p[0] = p[1] + p[2]

   def p_dist_el_ch1(self, p):
      'dist_el : opt_restrid LPAREN resrange_list RPAREN DASH LPAREN resrange_list RPAREN opt_chain LPAREN atom_and_weight RPAREN'
      p[0] = []
      atom = p[11]
      for res_id in p[3]:
         res_id.atom_name = atom[0]
         res_id.chain_id = p[9]
      for res_id in p[7]:
         res_id.atom_name = atom[1]
         res_id.chain_id = p[9]
      
      restr = DistanceRestraint(p[3],p[7], atom[2], atom[3],atom[4], p[1])
#      self.restrlist.append(restr)
      p[0].append(restr)


   def p_dist_el_ch2(self, p):
      'dist_el : opt_restrid LPAREN resrange_list RPAREN chain DASH LPAREN resrange_list RPAREN chain LPAREN atom_and_weight RPAREN'
      p[0] = []
      atom = p[12]
      for res_id in p[3]:
         res_id.atom_name = atom[0]
         res_id.chain_id = p[5]
      for res_id in p[8]:
         res_id.atom_name = atom[1]
         res_id.chain_id = p[10]
         
      restr = DistanceRestraint(p[3], p[8], atom[2], atom[3], atom[4], p[1])
      #print "residues!!!", p[3], p[8], p[11]
#      self.restrlist.append(restr)
      p[0].append(restr)
      
   
   def p_resrange_list_n(self, p):
      'resrange_list : resrange_el COMMA resrange_list'
      p[0] = p[3] + p[1]
       
   def p_resrange_list_1(self, p):
      'resrange_list : resrange_el'
      p[0] = p[1]  
         
   def p_oneres_el(self, p):
      'resrange_el : res'
      p[0] = [PDBId(p[1][0], p[1][1])]
      
   
   def p_resrange_el(self, p):
      'resrange_el : res DASH res'   
      p[0] = [PDBrange(p[1][0], p[1][1], p[3][0], p[3][1])]
       

######## atom distance specification with optional weight

   def p_atom_and_weight_0(self, p):
      'atom_and_weight : comp number opt_weight'
      at_el = ("", "", p[1], p[2], p[3])
      p[0] = at_el
      
   def p_atom_and_weight_1(self, p):
      'atom_and_weight : atom comp number opt_weight'
      at_el = (p[1], p[1], p[2], p[3], p[4])
      p[0] = at_el
      
   def p_atom_and_weight_2(self, p):
      'atom_and_weight : atom DASH atom comp number opt_weight'
      at_el = (p[1], p[3], p[4], p[5], p[6])
      p[0] = at_el


############SS Composition restraints 

   def p_ss_composition(self, p):
      'ss_composition : SS_COMPOSITION LPAREN ss_composition_list RPAREN'
      p[0] = p[3]

   def p_ss_composition_list_1(self, p):
      'ss_composition_list : opt_restrid INT COMP_LT ss_type COMP_LT INT opt_weight ss_composition_next'
      restraint = SSPercentageRestraint(p[4], p[2], p[6], weight=p[7], id=p[1])
      p[0] = [restraint]+p[8]

   def p_ss_type(self, p):
      '''ss_type : HELIX
                 | STRAND '''
      if   p[1].upper()=='HELIX':
        p[0] = ['H', 'G']
      elif p[1].upper()=='STRAND':
        p[0] = ['E', 'B']
      else: assert(False)

   def p_ss_composition_next_0(self, p):
      'ss_composition_next :'
      p[0] = []

   def p_ss_composition_next_1(self, p):
      'ss_composition_next : ss_composition_list'
      p[0] = p[1]
      
############Accessibility to surface - added by jmk
   def p_surface_acccess(self, p):
      'surface_access : SA LPAREN surface_acccess_list RPAREN'
      p[0] = p[3]

   def p_surface_acccess_list_0(self, p):
      'surface_acccess_list : '
      p[0] = []       
  
   def p_surface_acccess_list_1(self, p):
      'surface_acccess_list : surface_acccess_el surface_acccess_list'
      p[0] = p[1] + p[2]

   def p_surface_acccess_el(self, p):
      'surface_acccess_el : opt_restrid LPAREN resrange_list RPAREN opt_chain LPAREN atom_and_weight RPAREN'
      p[0] = []
      atom = p[7]
      for res_id in p[3]:
         res_id.atom_name = atom[0]
         res_id.chain_id = p[5]
      #print "residues!!!", p[3], "relation ", atom[2], "distval ", atom[3], "waga", atom[4], "p1", p[1]
      restr = SurfaceAccess(p[3], atom[2], atom[3],atom[4], p[1])
      #print "RR", restr
#      self.restrlist.append(restr)
      p[0].append(restr)
       


#############

############PointDistance - added by jmk
   def p_pointdist_acccess(self, p):
      'pointdist : PD LPAREN pointdist_list RPAREN'
      p[0] = p[3]

   def p_pointdist_list_0(self, p):
      'pointdist_list : '
      p[0] = []       
  
   def p_pointdist_list_1(self, p):
      'pointdist_list : pointdist_el pointdist_list'
      p[0] = p[1] + p[2]

   def p_pointdist_el(self, p):
      'pointdist_el : opt_restrid LPAREN resrange_list RPAREN opt_chain DASH LPAREN point RPAREN LPAREN atom_and_weight RPAREN'
      p[0] = []
      atom = p[11]
      for res_id in p[3]:
         res_id.atom_name = atom[0]
         res_id.chain_id = p[5]
      #print "residues!!!", p[3], "relation ", atom[2], "distval ", atom[3], "waga", atom[4]
      restr = PointDistance(p[3], p[8], atom[2], atom[3],atom[4], p[1])
      #print "RR", restr
#      self.restrlist.append(restr)
      p[0].append(restr)
      
   def p_point(self, p):
      'point : FLOAT COMMA FLOAT COMMA FLOAT'
      p[0] = [p[1], p[3], p[5]]
       


#############
   


############Symmetry distances - added by jmk
   def p_symmetry(self, p):
      'symmetry : SYM LPAREN symmetry_list RPAREN'
      p[0] = p[3]

   def p_symmetry_list_0(self, p):
      'symmetry_list : '
      p[0] = []       
  
   def p_symmetry_list_1(self, p):
      'symmetry_list : symmetry_el symmetry_list'
      p[0] = p[1] + p[2]

   def p_symmetry_el_ch1(self, p):
      'symmetry_el : opt_restrid LPAREN resrange_list RPAREN DASH LPAREN resrange_list RPAREN opt_chain '
      p[0] = []
      atom = p[11]
      for res_id in p[3]:
         res_id.atom_name = atom[0]
         res_id.chain_id = p[9]
      for res_id in p[7]:
         res_id.atom_name = atom[1]
         res_id.chain_id = p[9]
      
      restr = DistanceRestraint(p[3],p[7], atom[2], atom[3],atom[4], p[1])
      self.symrestrlist.append(restr)
      p[0].append(restr)
       
   def p_symmetry_el_ch2(self, p):
      'symmetry_el : opt_restrid LPAREN resrange_list RPAREN chain DASH LPAREN resrange_list RPAREN chain '
      p[0] = []
      atom = p[12]
      for res_id in p[3]:
         res_id.atom_name = atom[0]
         res_id.chain_id = p[5]
      for res_id in p[8]:
         res_id.atom_name = atom[1]
         res_id.chain_id = p[10]
         
      restr = DistanceRestraint(p[3], p[8], atom[2], atom[3], atom[4], p[1])
      #print "residues!!!", p[3], p[8], p[11]
      self.symrestrlist.append(restr)
      p[0].append(restr)
      
   
   def p_symresrange_list_n(self, p):
      'symresrange_list : resrange_el COMMA resrange_list'
      p[0] = p[3] + p[1]
       
   def p_symresrange_list_1(self, p):
      'symresrange_list : resrange_el'
      p[0] = p[1]  
         
   def p_symoneres_el(self, p):
      'symresrange_el : res'
      p[0] = [PDBId(p[1][0], p[1][1])]
      
   
   def p_resrange_el(self, p):
      'symresrange_el : res DASH res'   
      p[0] = [PDBrange(p[1][0], p[1][1], p[3][0], p[3][1])]
           



#############


############Accessibility restraints
   def p_access(self, p):
      'access : ACCESS LPAREN access_list RPAREN'
      p[0] = p[3]

   def p_access_list_0(self, p):
      'access_list : '
      p[0] = []

   def p_access_list_1(self, p):
      'access_list : access_el access_list'
      p[0] = p[1]+p[2]
      
   def p_access_el(self, p):
      'access_el : opt_restrid res opt_chain number opt_weight '
      res_id = PDBId(p[2][0],p[2][1],p[3])
      restr = AccessRestraint(res_id, p[4], p[4], p[5],p[1])
#      self.restrlist.append(restr)
      p[0] = [restr]

   def p_access_el_lr(self, p):
      'access_el : opt_restrid res opt_chain number DASH number opt_weight'
      res_id = PDBId(p[2][0],p[2][1],p[3])
      #print "#######", res_id
      restr = AccessRestraint(res_id, p[4],p[6], p[7],p[1])      
      #self.restrlist.append(restr)
      p[0] = [restr]

############# chains, atoms, residues, ids, etc.
   def p_opt_chain_empty(self, p):
      'opt_chain : '
      p[0] = ''

   def p_opt_chain(self, p): # sometimes we _can_ have chainID, but sometimes we _must_
      'opt_chain : chain'
      p[0] = p[1]

   def p_chain_id(self, p):
      'chain : CHAIN_ID'
      p[0] = p[1]

   def p_atom_id(self, p):
      'atom : ID'
      p[0] = p[1]
    
   def p_res_int(self, p):
      'res : INT'
      p[0] = ('',p[1])

   def p_res_id(self, p):
      'res : ID'
      p[0] = (p[1],None)
		
   def p_res_id_int(self, p):
      'res : ID INT'
      p[0] = (p[1],p[2])
		
#   def p_res_any(self,p):
#      'res : SIGN_ANY'
#      p[0] = (p[1],'')


   def p_opt_restrid_empty(self,p):
      'opt_restrid :'
      p[0] = ''

   def p_opt_restrid(self,p):
      'opt_restrid : RID'
      p[0] = p[1]

 

   def p_comp(self, p):
      ''' comp : COMP_LT
               | COMP_LE
               | COMP_EQ
               | COMP_GE
               | COMP_GT
      '''
      p[0] = p[1]


   def p_number_int(self, p):
      'number : INT'
      p[0] = float(p[1])
         
   def p_number_float(self, p):
      'number : FLOAT'
      p[0] = p[1]

    
   def p_error(self, p):
      #raise FormatError("Syntax error '%s' in line %d" % (p.value, p.lineno))
      print FormatError("Syntax error '%s' in line %d" % (p.value, p.lineno))


#########
 
      
if __name__ == '__main__':
   rp = RestraintsParser()
   (names,restlist) = rp.run(sys.argv[1])
   #(names,restlist) = rp.run('wiezy_dist_only.txt')
   for r in restlist:
      print r.name()
   
   print '*******names*******'
   for n in names:
      print n,names[n]

   print '*******IDs*********'
   for r in restlist:
      print r.get_id(), "\t" , r.reslist1
      if hasattr(r,'reslist2'): print r.reslist2
      print "sss", r.type, r.weight, r.check_residues
