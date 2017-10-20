#!/usr/bin/python
# Copyright 2006 Michal J. Gajda
# see license.txt

import sys, operator, re, itertools
from Bio.PDB.PDBParser import PDBParser
#from genesilico.structure.BrutePDBParser import BrutePDBParser as PDBParser
#import genesilico.structure.BrutePDBParser 
#########################
import tarfile
########################

def open_any_file(filename, mode="r"):
  """Open a file."""
  if filename == "-":
    if mode[0]=="r":
      return sys.stdin
    else: # Assuming "w"rite or "a"ppend
      return sys.stdout
  elif filename[-3:].lower() == ".gz":
    import gzip
    return gzip.open(filename, mode)
  elif filename[-3:].lower() == 'tar':
    return tarfile.open(filename, mode)
  elif filename[-3:].lower() == 'tar.gz':
    return tarfile.open(filename, "r:gz")
  else:
    return open(filename, mode)
  
class ModelSet(object):
# Should define just __iter__()
  pass
  

class CompositeModelSet(ModelSet):
  def __init__(self, modelsets):
    self.modelsets = modelsets
  def __iter__(self):
    return itertools.chain(*self.modelsets)
  def __len__(self):
    return reduce(operator.add,
                  map(len, self.modelsets),
                  0)

class PDBFilenameList(ModelSet):
  def __init__(self, pdblist):
    self.files = pdblist
  def __len__(self):
    return len(self.files)
  def __iter__(self):
    for fname in self.files:
      try:
        print fname
        file = open_any_file(fname, "r")
#        print "PATH:", sys.path, genesilico.__path__
#structure.BrutePDBParser.__path__
        parser = PDBParser(PERMISSIVE=False, QUIET=True)
        tmp_filename = fname
        if type(tmp_filename) != type(""):
          tmp_filename = tmp_filename.name
        structure = parser.get_structure(tmp_filename, file)
        file.close()
        structure.name     = tmp_filename
        structure.filename = tmp_filename
        yield structure
      except Exception, e:
        import traceback
        traceback.print_exc()
        sys.stderr.write("File %s doesn't exist or has wrong format.\n%s"
                           % (str(fname), str(e)))


class PDBDirfile(PDBFilenameList):
  def __init__(self, dirfilename, filesnames):
    self.dirfilename = dirfilename
    "zwraca liste nazw plikow wczytanych z pliku dirfile"
    
    #file = open_any_file(dirfilename, 'r') 
    #print "####", file
    #text = file.read()

    #file.close()
 
    #tmp_files = re.split('[ \t\n]+', text)
    #self.files = []
    #for file in tmp_files:
    #  print "^^^^^^^", file
    #  if file!='': self.files.append(file)
    #if len(self.files)>=1 and self.files[-1]=='':
    #  self.files=self.files[0:-1]
    self.files = filesnames


