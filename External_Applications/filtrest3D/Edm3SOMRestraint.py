# Copyright '2006 Marta Kaczor, Michal J. Gajda, Anastasia Bakulina
# See license.txt

import PModel  #getting access to stride_dict
from Restraint import Restraint
import sys
import tempfile
import os
from mapconverter import convert

class Edm3SOMRestraint(Restraint):
    'docking to an electron density map, using 3SOM program'
    
    type=["EDM","3SOM"]
    
    def __init__(self, edm_filename, id):
       #self.edm_filename = edm_filename
       if edm_filename.split('.')[-1].lower() == 'ccp4':
          self.edm_filename = edm_filename
       else: #convertion is needed
          try: (fd, converted_filename) = tempfile.mkstemp()     
          except Exception, e:
             print >> sys.stderr, "*** tempfile.mkstemp() failed: "+str(e)+" ***"
             self.edm_filename = ''
             return 0
          os.close(fd)   
          if convert(edm_filename, 'ccp4', converted_filename) != 0: 
             print >> sys.stderr, '*** Cannot convert map file!'
             self.edm_filename = ''
          else: self.edm_filename = converted_filename
       self.voxel_size = 3
       self.resolution = 20
       self.id = id
       
    def check(self, pmodel):
       if self.edm_filename == '':
          print >> sys.stderr,'*** Error in EDM restraint!'
          return 0
       try: (fd, out_3som) = tempfile.mkstemp()
       except Exception, e:
          print >> sys.stderr, "*** tempfile.mkstemp() failed: "+str(e)+" ***"
          return 0
       os.close(fd)
       callstring = '~/bin/3som/3som_silent '+self.edm_filename+' '+pmodel.filename+' '+str(self.voxel_size)+' '+str(self.resolution)+' '+out_3som
       ret = os.system(callstring)                                                                           
       if ret != 0: 
          print >> sys.stderr, '*** error in 3som program!' 
          return 0

       file = open(out_3som, 'r')
       line = file.readline()
       items = filter(None,line.split(' '))
       percents = items[5][0:-2].strip() #remove ',' and spaces
       score = (1 - float(percents))*self.PENALTY_MAX
       return score
       
    def autoname (self):
       return self.edm_filename+' '+str(self.voxel_size)+' '+str(self.resolution)
       
       
