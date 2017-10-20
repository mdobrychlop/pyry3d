# wrapper to Situs/conformat
# Copyright 2005 Anastasia Bakulina
# see license.txt
# usage: mapconverter.py filename type
# input:                ccp4 situs pdb
# output(def. by type): ccp4 situs
import sys
import os
import tempfile

def convert(filename, out_type, out_file = ''):
   supported_input  = ['ccp4','situs','pdb']
   supported_output = ['ccp4','situs']
   conformat = '~/bin/Situs/src/conformat_silent'
   pdblur = '~/bin/Situs/src/pdblur_silent'
   #check input
   try:
      testopen = open(filename)
   except IOError:
      print 'Input file not found!'
      return 1   
   out_type = out_type.lower()
   if out_type not in supported_output:
       print 'unsupported format for output. Now support only:',
       for f in supported_output: print f,
       print '\n'
       return 1
       
   fn = filename.split('.')
   
   if len(fn) < 2:
      print 'input file should has an extension'
      return 1
      
   dot = '.'
   in_name = dot.join(fn[0:-1])
   in_ext = fn[-1].lower()
   if in_ext not in supported_input:
      print 'unsupported format for input. Now support only:',
      for f in supported_input: print f,
      print '\n'
      return 1
   
   if in_ext == out_type:
      print 'I have nothing to do'
      return 0
   
   if out_file != '': outfilename = out_file
   else:              outfilename = in_name + '.' + out_type   
   space = ' '
   output = '> /dev/null'
   if in_ext == 'ccp4' and out_type == 'situs':
      callstring = space.join((conformat,filename,outfilename,str(3),output))
      ret = os.system(callstring)
      if ret != 0:
         print 'error in conformat!'
         return 1
   elif in_ext == 'situs' and out_type == 'ccp4':
      callstring = space.join((conformat,filename,outfilename,str(8),output))
      ret = os.system(callstring)
      if ret != 0:
         print 'error in conformat!'
         return 1
   elif in_ext == 'pdb' and out_type == 'situs':
      callstring = space.join((pdblur, filename, outfilename,str(2),output))#let voxel size = 2
      ret = os.system(callstring)
      if ret != 0:
         print 'error in pdblur!'
         return 1
   elif in_ext == 'pdb' and out_type == 'ccp4':
      try: (fd, tmpfilename) = tempfile.mkstemp()
      except Exception, e:
         print >> sys.stderr, "*** tempfile.mkstemp() failed: "+str(e)+" ***"         
      os.close(fd)
      callstring = space.join((pdblur, filename, tmpfilename,str(2),output))#let voxel size = 2   
      ret = os.system(callstring)
      if ret != 0:
         print 'error in pdblur!'
         return 1         
      callstring = space.join((conformat,tmpfilename,outfilename,str(8),output))      
      ret = os.system(callstring)
      if ret != 0:
         print 'error in conformat!'
         return 1
      os.remove(tmpfilename)
   return 0

if __name__ == '__main__':
   if len(sys.argv) != 3: 
      print 'Usage: mapconverter.py input_file output_type'
   else:
      if convert(sys.argv[1], sys.argv[2]) == 0:
         print 'Done.'
      else:
         print 'error.'   
