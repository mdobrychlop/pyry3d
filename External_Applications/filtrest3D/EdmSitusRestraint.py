# Copyright '2006 Marta Kaczor, Michal J. Gajda, Anastasia Bakulina
# See license.txt

import PModel  #getting access to stride_dict
from Restraint import Restraint
import sys
import tempfile
import os
from mapconverter import convert

class EdmSitusRestraint(Restraint):
    'docking to an electron density map, using simplified Situs program'
    
    type=["EDM","Situs"]
    
    def __init__(self, edm_filename, id): #maybe I'll need more parameters
       if edm_filename.split('.')[-1].lower() == 'situs':
          self.edm_filename = edm_filename
       else: #convertion is needed
          try: (fd, converted_filename) = tempfile.mkstemp()
          except Exception, e:
             print >> sys.stderr, "*** tempfile.mkstemp() failed: "+str(e)+" ***"
             self.edm_filename = ''
             return 0
          os.close(fd)
          if convert(edm_filename, 'situs', converted_filename) != 0:
             print >> sys.stderr, '*** Cannot convert map file!'
             self.edm_filename = ''
          else: self.edm_filename = converted_filename
       self.id = id
       
    def check(self, pmodel):
       if self.edm_filename == '':
          print >> sys.stderr,'*** Error in EDM restraint!'
          return 0
    
       try: (fd1, out_situs) = tempfile.mkstemp()
       except Exception, e:
          print >> sys.stderr, "*** tempfile.mkstemp() failed: "+str(e)+" ***"
          return 0
       os.close(fd1)
       try: (fd2, err_situs) = tempfile.mkstemp()
       except Exception, e:
          print >> sys.stderr, "*** tempfile.mkstemp() failed: "+str(e)+" ***"
          return 0
       os.close(fd2)

       callstring = '~/bin/Situs/qrange_silent '+self.edm_filename+' '+pmodel.filename+ ' 1>' + out_situs + ' 2>' + err_situs
       ret = os.system(callstring)                                                                           
       if ret != 0:
          print >> sys.stderr, "*** Error in Situs program! " # I don't now what ret or ret>>8 means, so I don't output it
          return 0 
       file = open(out_situs, 'r')
       line = file.readline()
       items = filter(None,line.split(' '))
       while items[0] != 'qrange>':  
          line = file.readline()          
          items = filter(None,line.split(' '))
#          while items == []: 
#             line = file.readline()          
#             items = filter(None,line.split(' '))
       return (1-float(items[2]))*self.PENALTY_MAX # or items[1]?
       
    def autoname (self):
       return self.edm_filename
       
       
