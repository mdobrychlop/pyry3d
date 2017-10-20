#!/usr/bin/env python
from __future__ import division
# -*- coding: utf-8 -*-
#
# 
# www.genesilico.pl 
#

#creates ranked 3D models of macromoleular complexes 
#based on experimental restraints and a whole complex shape.


__author__ = "Joanna M. Kasprzak"
__copyright__ = "Copyright 2010, The PyRy3D Project"
__credits__ = ["Janusz Bujnicki"]
__license__ = "GPL"
__version__ = "0.1.0"
__maintainer__ = "Joanna Kasprzak"
__email__ = "jkasp@amu.edu.pl"
__status__ = "Prototype"


import sys

class Outfolder:
    """
    Represents log data returned by the program both to view on screen and
    write to file
    """
    def __init__(self):
        """
            filename is optional
        """
        self.outdirname = ''
        
    def set_dirname(self, dirname):
        """
        """
        self.outdirname = dirname


class Log:
    """
    Represents log data returned by the program both to view on screen and
    write to file
    """
    def __init__(self):
        """
            filename is optional
        """
        self.contents = [] ## list
        #self.__add_header()
        self.write_to_stderr = True
        self.print_all = False
        self.file_name = ""

    def __del__(self):
        """
            It is called before an object of this class is destroyed
        """
        if self.contents:
            self.write_file()

    def set_filename(self, name):
        """
            Sets the name of the log file
        
        Parameters:
        ----------
            name    : name of log file if different than pyry.log
        """
        self.file_name = name

    def write_message(self, message):
        """
            Adds a string message to the log
        
        Parameters:
        -----------
            message :   info for program user
        """
        self.contents.append(message+'\n')
        if self.print_all:
            print message
            
    def write_error(self, error):
        """
            Adds an Exception to the log
        
        Parameters:
        -----------
            error   : ERROR message
        """
        self.contents.append('\nERROR: '+str(error)+'\n')
        # also writing errors to stderr, so they are visible
        # in the console and shell
        if self.write_to_stderr:
            sys.stderr.write(str(error)+'\n')
    
    def write_file(self, text=False):
        """
            Writes all log messages to a file
        """
        f = open(self.file_name,'a')
        if text: f.writelines(text)
        else: f.writelines(self.contents)
        f.close()
        self.__clear_file()

    def __clear_file(self):
        """
            Clears all log messages
        """
        self.contents = []
    
    def __add_header(self):
        """
            Generates the header of log file
        """
        header = """\
********************************************************************

    PyRy3D 
    program for macromolecular docking into cryoEM maps
    
    (c) 2010 by Joanna M. Kasprzak, Janusz M. Bujnicki
    
    usage: python pyry3d.py --help
    
    maintainer: jkasp@amu.edu.pl
    
*********************************************************************   

        \n"""
        self.contents.append(header)
        
        
logfile = Log()
trafl = Log()
movehist = Log()
outfolder = Outfolder()
