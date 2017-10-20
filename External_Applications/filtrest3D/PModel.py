# -----------------------------------------------------------------------------
# PModel.py
#
# 
# -----------------------------------------------------------------------------
# Copyright 2003-2006 Marta Kaczor, Michal J. Gajda, Anastasia Bakulina
# see license.txt


import os
from stat import *
import re
#from Bio.PDB import *

#####
from Bio.PDB import PDBParser
#from Bio.SCOP.Raf import to_one_letter_code
from Modules.Constans.Constans  import to_one_letter_code

#####

from Bio.PDB.StructureBuilder import StructureBuilder
#from StructureBuilder import StructureBuilder
from string import split
#from Numeric import array, Float0

from PDBId import PDBId

import math
import tempfile
import sys
import string

def pick_binary(binaries_list):
  # Find first path from a list that exists and return it.
  # Does full shell expansion, but no globbing.
  for binary_path in binaries_list:
    binary_path = os.path.expandvars(os.path.expanduser(binary_path))
    if os.path.exists(binary_path):
#      sys.stderr.write("Picked binary %s\n" % binary_path)
      return binary_path
  
stride_path=pick_binary(["~/bin/stride",
                         "./stride",
                         "/usr/local/bin/stride",
                         "/usr/bin/stride"])

dssp_path=pick_binary(["dssp",
                       "~/bin/dssp",
                       "/usr/local/bin/dssp",
                       "/usr/bin/dssp"])
sable_path="run_sable.sh"

DIST_NO_RESIDUE_PENALTY = 100.0
KNOT_PENALTY = 100.0

ACCESS_NO_RESIDUE_PENALTY = 100.0

SABLE_FAILED_TO_RUN_PENALTY = 1000.0

MODELS_IS_SAME = 1

global f_punish
f_punish = 0

Ala_X_Ala = {
   'A' : 110.2, 'D' : 144.1, 'C' : 140.4, 'E' : 174.7, 'F' : 200.7, 
   'G' :  78.7, 'H' : 181.9, 'I' : 185.0, 'K' : 205.7, 'L' : 183.1, 
   'M' : 200.1, 'N' : 146.4, 'P' : 141.9, 'Q' : 178.6, 'R' : 229.0, 
   'S' : 117.2, 'T' : 138.7, 'V' : 153.7, 'W' : 240.5, 'Y' : 213.7, 
   }

MAXIMAL_ACC = {
   'A' : 106, 'C' : 135, 'D' : 163, 'E' : 194, 'F' : 197,
   'G' :  84, 'H' : 184, 'I' : 169, 'K' : 205, 'L' : 164,
   'M' : 188, 'N' : 157, 'P' : 136, 'Q' : 198, 'R' : 248,
   'S' : 130, 'T' : 142, 'V' : 142, 'W' : 227, 'Y' : 222
   }



mkstemp_prefix = 'filtrest3d.tmp.file'
mkstemp_suffix = '.tmp'

class EmptyStructure(Exception):
  def __init__(self, structure):
    self.structure = structure
  def __str__(self):
    return "Empty structure " + self.structure

class PModel(object):
   '''
   klasa reprezentujaca model bialka, sprawdzajaca spelnienie wiezow i pamietajaca wyniki tych obliczen

   do stworzenia obiektu wymaga nazwy pliku w formacie PDB
   '''


   def __init__(self, structure, use_sable=False):
      'param filename - nazwa pliku w formacie PDB'
      
      self.structure = structure
      try:
        self.model = structure[0] # structure[0], structure is returned by get_structure 
      except:
        raise EmptyStructure, structure.name
      self.filename = structure.name

      self.sable_acc_dict = {}
      self.sable_acc_res = []
      
      self.sable_query = None   #used
      self.sable_translate = None #used
      self.sable_res = None
      self.check_results = {} # dict of results of checking restraints
      self.knotted = None
      self.first_chain_id = ''#identifier of the first chain of the model, needs when chainID is omitted in a restraint
      self.multichain = 0

      if (use_sable):
         self.sable_query = ''
         self.sable_translate = {}
         sable_ind = 0

      if self.first_chain_id == ' ': self.first_chain_id = ''
      
      if len(self.model.get_list()) > 1: self.multichain = 1 
      
      for chain in self.model.get_list():
         for residue in chain.get_list():
            res_name = string.strip(residue.get_resname())
            res_num = residue.get_id()[1]
            if len(res_name)==3: 
               if res_name in to_one_letter_code:
                  res_name = to_one_letter_code[res_name]
            if use_sable and len(res_name)==1:
               self.sable_query += res_name
               self.sable_translate[sable_ind] = res_num
               sable_ind += 1
            #end if
         #end for
      #end for

   def delPDBStructures(self):
      del self.model
      return

   def get_iterator(self):
      return self.model.get_iterator()
           
   

   def stride(self):
      '''
      uruchamia stride, zapisuje wyniki (dict: res_num -> {'H','G','I','E','B','b','C','T'}) w self.stride_dict
      zwraca 0 - sukces, -1 - blad
      '''
      
      try: (fd, out_tmp1) = tempfile.mkstemp(prefix = mkstemp_prefix, suffix = mkstemp_suffix)
      except Exception, e:
         print >> sys.stderr, "*** tempfile.mkstemp() failed: "+str(e)+" ***"
         return -1
      os.close(fd)
      
      
      try: (fd, out_tmp2) = tempfile.mkstemp(prefix = mkstemp_prefix, suffix = mkstemp_suffix)
      except Exception, e:
         print >> sys.stderr, "*** tempfile.mkstemp() failed: "+str(e)+" ***"
         os.remove(out_tmp1)
         return -1
      os.close(fd)
      
      try: (fd, err_tmp) = tempfile.mkstemp(prefix = mkstemp_prefix, suffix = mkstemp_suffix)
      except Exception, e:
         print >> sys.stderr, "*** tempfile.mkstemp() failed: "+str(e)+" ***"
         os.remove(out_tmp1)
         os.remove(out_tmp2)
         return -1
      os.close(fd)
      
      polecenie = "%s %s 1>%s 2>%s" % (stride_path, self.filename, out_tmp1, err_tmp)
      # moj stride - zmianilam mu w funkcji die() exit(1) na exit(FAILURE) - zeby nie zwracal zawsze 1
      ret = os.system(polecenie)
      ret = ret >> 8
      
      myret = 0
      if ret != 1: #niestety, dla stride SUCCESS to 1, FAILURE to 0 :)
         print >> sys.stderr, "*** Error: failed to run stride on file "+self.filename+" ***"
         ferr = open(err_tmp)
         text = ferr.read()
         ferr.close()
         if len(text)>0:
            if text[-1]=='\n': text = text[0:-1]
            print >> sys.stderr, "*** *** stride.err: ["+text+"] *** ***"
         #end if
         myret = -1
      else:
         os.system("perl -pe 's/\x01@/\x0a/go' <%s >%s" % (out_tmp1, out_tmp2))
         # kolejna mila cecha stride :) - dla jednego pliku z moich testowych potrzebna byla taka zamiana
         # - ale tylko jak stride dzialal na drugs - ???
         file = open(out_tmp2, 'r')
         lines = file.readlines()
         file.close()
      
         regexp = re.compile('^ASG .*')
         for line in lines:
            if regexp.match(line): 
               resnum = int(line[11:15])
               chain_id = line[9]
               if chain_id == '-': chain_id = ''
               resid = PDBId('',resnum,chain_id)
               res = resid.find1res(self)
               res.secondary_structure = line[24]
         #end for
         myret = 0
      #end if
      
      os.remove(err_tmp)
      os.remove(out_tmp1)
      os.remove(out_tmp2)
      return myret
   
   
   def dssp(self):
      'uruchamia dssp, zapisuje wynik w self.dssp_acc_dict (dict: res_num -> (100 * ACC z DSSP(x) / Ala_X_Ala(x)))'
      'zwraca 0 -udalo sie, -1 - blad'

      try: (fd, out_tmp) = tempfile.mkstemp(prefix = mkstemp_prefix, suffix = mkstemp_suffix)
      except Exception, e:
         print >> sys.stderr, "*** tempfile.mkstemp() failed: "+str(e)+" ***"
         return -1
      os.close(fd)
      
      try: (fd, err_tmp) = tempfile.mkstemp(prefix = mkstemp_prefix, suffix = mkstemp_suffix)
      except Exception, e:
         print >> sys.stderr, "*** tempfile.mkstemp() failed: "+str(e)+" ***"
         os.remove(out_tmp)
         return -1
      os.close(fd)
      
      ### I havn't dssp installed on my system, so I'll just use out file:
      ###out_tmp = '/home/nastya/filtrest3d_new/dssp.out'
      polecenie = "%s %s 1>%s 2>%s" % (dssp_path, self.filename, out_tmp, err_tmp)
      ret = os.system(polecenie)
      ret = ret >> 8
      ret = 0
      
      myret = 0
      if ret != 0 or not os.path.exists(out_tmp):
         print >> sys.stderr, "*** Error: failed to run dssp on file "+self.filename+" ***"
         ferr = open(err_tmp)
         text = ferr.read()
         ferr.close()
         if len(text)>0:
            if text[-1]=='\n': text = text[0:-1]
            print >> sys.stderr, "*** *** dssp.err: ["+text+"] *** ***"
         #end if
         myret = -1
      else:
         myret = 0
         file = open(out_tmp, 'r')
         lines = file.readlines()
         file.close()

         regexp1 = re.compile('[ |\t]*#[ |\t]*RESIDUE.*')
      
         i = -1
         llen = len(lines)
         while (i < llen):
            i+=1
            line = lines[i]
            if not regexp1.match(line):
               continue
            b_acc = string.find(line, "ACC")
            e_acc = b_acc + 3
            b_aa = string.find(line, "AA")
            e_aa = b_aa + 2
            if b_aa <= 0 or b_acc <= 0:
               myret = -1
            break
         #end while
         
         if myret != -1:
            for j in range(i+1, llen):
               try: 
                 line=lines[j]
                 res_num = int(line[6:10])
                 chain_id = line[11]
                 if chain_id == ' ': chain_id = ''
                 acc = int(line[b_acc:e_acc])
                 aa_name = string.strip(line[b_aa:e_aa])

                 dssp_value = float(100*acc)/Ala_X_Ala[aa_name]
                 if dssp_value > 100.0:
                    dssp_value = 100.0
                 
                 
                 resid = PDBId('',res_num,chain_id)
                 res = resid.find1res(self)
                 res.dssp = dssp_value 
               
               except:
                 sys.stderr.write("Cannot parse DSSP output line:\n%s\n" % line)
            #end for
         #end if myret != -1
      #end if
      
      os.remove(err_tmp)
      os.remove(out_tmp)
      return myret



   def __run_sable_sh__(self, sable_tmpdir, outf, errf):
      'uruchamia sable '
      try: (fd, errf_sh) = tempfile.mkstemp(prefix = mkstemp_prefix, suffix = mkstemp_suffix, dir=sable_tmpdir)
      except Exception, e:
         print >> sys.stderr, "*** tempfile.mkstemp() failed: "+str(e)+" ***"
         return -1
      os.close(fd)

      #ret = os.system('./run.sable.sh %s ' % (self.sable_query))
      ret = os.system('%s %s %s %s %s 1>/dev/null 2>%s' % (sable_path, self.sable_query, outf, errf, sable_tmpdir, errf_sh))
      ret = ret >> 8
      comm = ''
      if ret != 0:
         try: f = open(errf)
         except:
            f = open(errf_sh)
            text = f.read()
            f.close()
            os.remove(errf_sh)
            if len(text)>0 and text[-1]=='\n': text = text[0:-1]
            if len(text)>0: comm = 'failed to run \'run.sable.sh\': '+text
            else: comm = 'failed to run \'run.sable.sh\''
            return -1

         text = f.read()
         f.close()
         if len(text)>0 and text[-1]=='\n': text = text[0:-1]
         comm = 'failed to run sable'
         if len(text)>0: comm += ': "%s"' % text
         print >> sys.stderr, '*** FAILED to run sable',
         if len(text)>0: print >> sys.stderr, ': "%s" ***' % text
         else: print >> sys.stderr, ' ***'
      os.remove(errf_sh)
      return comm



   def __make_sable_acc_dict__(self, sableout):
      
      f = open(sableout)
  
      lines = f.readlines()
      f.close()
   
      print lines

      comm = 'unexpected sable output file format'
      
      if len(lines)<10 : return comm + '(1)'
      i = 0
      count = len(lines)
      empty_line = re.compile('^( |/t|/n)*$')

      # (ewentualne) poczatkowe puste linie
      while i<count and empty_line.match(lines[i]): i+=1
      if i>=count: return comm + '(2)'

      # linia "Query: ..."
      if not re.match('^[Q|q][u|U][e|E][r|R][y|Y]:.*$', lines[1]): return comm + '(3)'
      i+=1
      
      # linia "SECTION_SA"
      if not re.match('^SECTION_SA.*$', lines[i]): return comm + '(4)'
      i+=1

      # (ewentualne) puste linie
      while i<count and empty_line.match(lines[i]): i+=1
      if i>=count: return comm + '(5)'

      # cztery linie objasniajace, co jest w parsowanym pliku ponizej
      i+=4

      # (ewentualne) puste linie
      while i<count and empty_line.match(lines[i]): i+=1
      if i>=count: return comm + '(6)'
      
      # linie "> nr_rez          nr_rez "
      #          nazwy reziduow
      #          wspolczynniki przewidywanej solwatacji
      #          wspolczynniki pewnosci powyzszych
      
      res_count = 0
      end_section = re.compile('^END_SECTION.*$')
      while i<count:
         if end_section.match(lines[i]): break
         if lines[i][0]!='>': return comm + '(6)'
         if i+3>=count: return comm + '(7)'
         line1=lines[i+1][0:-1]
         line2=lines[i+2][0:-1]
         line3=lines[i+3][0:-1]
         i+=4
         for j in range(0, len(line1)):
            if line1[j]==' ': continue
            if line1[j]!=self.sable_query[res_count]: return comm+ '(8)'
            #print >> sys.stderr, res_count,'->',line1[j], line2[j], line3[j]
            try:
               if len(line2) <= j or len(line3) <= j: return comm+ '(9)'
               res = PDBId("",self.sable_translate[res_count],"").find1res(self)
               res.sable = (int(line2[j]), int(line3[j]))
               print res.get_resname(), res.get_id()[1], res.sable
            except Exception, e:
               return comm+ 'something wrong in __make_sable_acc_dict__'
            res_count+=1
         #end for j in range(0, len(line1))
      #end while i<count
      if res_count != len(self.sable_query): return comm+ '(10)'
      return ''


   
   def checkAccess_sable(self):
      'uruchamia sable i ocenia wyniki solwatacji sable z solwatacja odczytana z dssp'
      'dziala dla nowej wersji sable (sable2)- tej z "confidence level"'
      'zwraca: 0 - udalo sie, -1 - blad'
      
      try:
         sable_tmpdir = tempfile.mkdtemp(prefix = mkstemp_prefix, suffix = mkstemp_suffix)
      except Exception, e:
         print >> sys.stderr, '*** tempfile.mkdtemp() failed: '+str(e)+' ***'
         return -1
      
      try: (fd, outf) = tempfile.mkstemp(prefix = mkstemp_prefix, suffix = mkstemp_suffix, dir = sable_tmpdir)
      except Exception, e:
         print >> sys.stderr, "*** tempfile.mkstemp() failed: "+str(e)+" ***"
         return -1
      os.close(fd)
      try: (fd, errf) = tempfile.mkstemp(prefix = mkstemp_prefix, suffix = mkstemp_suffix, dir = sable_tmpdir)
      except Exception, e:
         print >> sys.stderr, "*** tempfile.mkstemp() failed: "+str(e)+" ***"
         return -1
      os.close(fd)
      

###      comm = self.__run_sable_sh__(sable_tmpdir, outf, errf)

###      if comm != '':
###         self.sable_res = (SABLE_FAILED_TO_RUN_PENALTY, comm)
###         os.system('rm -fr %s' % sable_tmpdir)
###         return 0
      
      comm = self.__make_sable_acc_dict__(sable_tmpdir+'/'+'OUT_SABLE_RES')
      print comm
      if comm != '':
         self.sable_res = (SABLE_FAILED_TO_RUN_PENALTY, comm)
         os.system('rm -fr %s' % sable_tmpdir)
         return 0

      os.system('rm -fr %s' % sable_tmpdir)
      return 0


   def delDsspStructures(self):
      del self.sable_acc_dict
      del self.sable_query
      del self.sable_translate
   
   
   def summarize(self): # we suppose it will be called only once
      sumtoken = 'sum'
      if self.check_results == {}: return 0 # empty
      sumresults = {(sumtoken,):0}

      for checkresult in self.check_results:
         sumname = (sumtoken,)
         sumresults[sumname] += self.check_results[checkresult] # add total sum
            
#         for restrtype in checkresult[0]:  # for example, we sum all "Distance",then all "Distance" "Fuzzy"
#            sumname += (restrtype,)
#            if sumresults.get(sumname,None) == None: #new type of restr 
#               sumresults.update({sumname:0.0})
#            sumresults[sumname] += self.check_results[checkresult]
   
      self.check_results.update(sumresults)


#end class PModel

#
   def _parse_coordinates(self, coords_trailer):
      "Parse the atomic data in the PDB file."
      local_line_counter=0
      structure_builder=self.structure_builder
      current_model_id=0
      current_chain_id=None
      current_segid=None
      current_residue_id=None
      current_resname=None
      structure_builder.init_model(current_model_id)
      for i in range(0, len(coords_trailer)):
          line=coords_trailer[i]
          record_type=line[0:6]
          global_line_counter=self.line_counter+local_line_counter+1
          structure_builder.set_line_counter(global_line_counter)
          if(record_type=='ATOM  ' or record_type=='HETATM'):
              fullname=line[12:16]
              # get rid of whitespace in atom names
              split_list=split(fullname)
              if len(split_list)!=1:
                  # atom name has internal spaces, e.g. " N B ", so
                  # we do not strip spaces
                  name=fullname
              else:
                  # atom name is like " CA ", so we can strip spaces
                  name=split_list[0]
              altloc=line[16:17]
              resname=line[17:20]
              chainid=line[21:22]
              resseq=int(split(line[22:26])[0])   # sequence identifier   
              icode=line[26:27]           # insertion code
              if record_type=='HETATM':       # hetero atom flag
                  if resname=="HOH" or resname=="WAT":
                      hetero_flag="W"
                  else:
                      hetero_flag="H"
              else:
                  hetero_flag=" "
              residue_id=(hetero_flag, resseq, icode)
              # atomic coordinates
              x=float(line[30:38]) 
              y=float(line[38:46]) 
              z=float(line[46:54])
              coord=array((x, y, z), Float0)
              # occupancy & B factor
              
              # we don't need it:
              try:
                occupancy=float(line[54:60])
              except ValueError:
                occupancy=1.0
              try:
                bfactor=float(line[60:66])
              except ValueError:
                bfactor=1.0
              #######
              
              segid=line[72:76]
              
              
              if current_segid!=segid:
                  current_segid=segid
                  structure_builder.init_seg(current_segid)
              if current_chain_id!=chainid:
                  current_chain_id=chainid
                  structure_builder.init_chain(current_chain_id)
                  current_residue_id=residue_id
                  current_resname=resname
                  try:
                      structure_builder.init_residue(resname, hetero_flag, resseq, icode)
                  except PDBConstructionException, message:
                      self._handle_PDB_exception(message, global_line_counter)
              elif current_residue_id!=residue_id or current_resname!=resname:
                  current_residue_id=residue_id
                  current_resname=resname
                  try:
                      structure_builder.init_residue(resname, hetero_flag, resseq, icode)
                  except PDBConstructionException, message:
                      self._handle_PDB_exception(message, global_line_counter) 
              # init atom
              try:
                  structure_builder.init_atom(name, coord, bfactor, occupancy, altloc, fullname)
              except PDBConstructionException, message:
                  self._handle_PDB_exception(message, global_line_counter)
          elif(record_type=='ANISOU'):
              anisou=map(float, (line[28:35], line[35:42], line[43:49], line[49:56], line[56:63], line[63:70]))
              # U's are scaled by 10^4 
              anisou_array=(array(anisou, Float0)/10000.0).astype(Float0)
              structure_builder.set_anisou(anisou_array)
          elif(record_type=='ENDMDL'):
              current_model_id=current_model_id+1
              structure_builder.init_model(current_model_id)
              current_chain_id=None
              current_residue_id=None
          elif(record_type=='END   ' or record_type=='CONECT'):
              # End of atomic data, return the trailer
              self.line_counter=self.line_counter+local_line_counter
              return coords_trailer[local_line_counter:]
          elif(record_type=='SIGUIJ'):
              # standard deviation of anisotropic B factor
              siguij=map(float, (line[28:35], line[35:42], line[42:49], line[49:56], line[56:63], line[63:70]))
              # U sigma's are scaled by 10^4
              siguij_array=(array(siguij, Float0)/10000.0).astype(Float0)   
              structure_builder.set_siguij(siguij_array)
          elif(record_type=='SIGATM'):
              # standard deviation of atomic positions
              sigatm=map(float, (line[30:38], line[38:45], line[46:54], line[54:60], line[60:66]))
              sigatm_array=array(sigatm, Float0)
              structure_builder.set_sigatm(sigatm_array)
          local_line_counter=local_line_counter+1
      # EOF (does not end in END or CONECT)
      self.line_counter=self.line_counter+local_line_counter
      return []

 
