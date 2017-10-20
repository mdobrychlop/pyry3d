#!/usr/bin/env python
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

import sys, os #, traceback
from optparse                                          import OptionParser
from External_Applications.filtrest3D.RestraintsParser import *
from External_Applications.filtrest3D.PModel           import PModel, EmptyStructure
from External_Applications.filtrest3D.ModelSet         import *
from DataDownloader                                    import * 
from Modules.Error.Errors                              import InputError
from Modules.Input.InputGenerator                      import *
                   

def get_args(inp):
    """
        uses OptionParser to get command line arguments: sequences, structures,
        restraints and a density map file    
    Parameters:
    -----------
        inp             : Input object
    Returns:
    --------
        sequences       : user defined seqs as Bio objects
        structures      : user defined structures as Bio objects
                          will be transformed into Complex_component instance
        restraints      : user defined restraints, will become
                          Distance_Interactions instances
        complex_map     : complex map data which will be transformed
                          into Complex_component instance
        
        """
    parser = OptionParser()
    parser.add_option("-s", "--seq_file", dest="seq_filename",
                  help="input file with sequences in FASTA format",\
                       metavar="FILE")
    parser.add_option("-d", "--dirfile", dest="dirname",
                  help="input file with structure", metavar="FILE")
    parser.add_option("-r", "--restr_file", dest="restr_filename",
                  help="input file with restraints in Filtrest3D format", \
                       metavar="FILE")
    parser.add_option("-m", "--map_file", dest="map_filename",
                  help="input file with complex map", metavar="FILE")
    parser.add_option("-x", "--saxs_file", dest="saxs_filename",
                  help="input file with complex abinitio reconstruction", metavar="FILE")
    parser.add_option("-a", "--auto", dest="auto",
                  help="generate input files automatically from structures", metavar="FILE")
    parser.add_option("-c", "--config", dest="config",
                  help="config file with simulation parameters")
    parser.add_option("-o", "--output", dest="output",
                  help="write output to file")
    parser.add_option("-t", "--traflfile", dest="traflfile",
                  help="write output to trajectory file")
    parser.add_option("-y", "--.dat file with experimental curve from SAXS", dest="curve_filename",
                  help="Verify discrepancy between theoretical and experimental curves from SAXS/SANS. ")
    parser.add_option("-f", "--to full atom model", dest="fullatom",
                  help="write best model in fullatom representation")
    parser.add_option("-v", "--save history of moves", dest="movehistory",
                  help="save all moves into a file")
    parser.add_option("-e", "--save plot with complexes scores", dest="score_plot",
                  help="save plot with complexes scores")
    parser.add_option("--fast", action="store_true", dest="fast", default=False,
                  help="run optimized version, option -v is then disabled")
    #other options???
    (opts, args) = parser.parse_args()
    modelsets = []


#------------optimization-----------------------------------------------
    inp.optimized = opts.fast

#---------------get sequences-------------------------------------
    if opts.seq_filename != None:
        seqs = Sequences()
        inp.sequences = seqs.get_seqs_data(opts.seq_filename)

#---------------get structures-------------------------------------
    if opts.dirname != None:
        structs = Structures()
        inp.structures = structs.get_structures(opts.dirname)
        modelsets.append(PDBDirfile(opts.dirname, structs.pdb_files))
        modelset=CompositeModelSet(modelsets)
    if opts.dirname == None:   
        inp.structures = []
        #raise InputError("No structure folder given, I have nothing to do!")
        
#---------------get restraints-------------------------------------
    if opts.restr_filename != None:
        restrs = PyryRestraints()
        inp.restraints = restrs.get_restraints(opts.restr_filename) #, modelset)
        if opts.restr_filename and inp.restraints == []:
            if inp.symrestraints == []:
                raise InputError("Format of file with restraints is wrong")
    else:
        restraints = []       
    if opts.restr_filename == None:
        #if restraint option is not used, than use empty list of restraints!
        inp.restraints = []
        
#---------------get density map or saxs shape----------------------------------
    if opts.map_filename != None:   
        map = ShapeDescriptor()
        inp.mapfile = map.get_density_map(opts.map_filename)
    if opts.saxs_filename != None:
        saxs = ShapeDescriptor()
        inp.saxsfile = saxs.get_saxs_shape(opts.saxs_filename)
    if opts.curve_filename != None:   
        map = ShapeDescriptor()
        inp.curvefile = map.get_curve(opts.curve_filename)
    elif (opts.map_filename == None) and (opts.saxs_filename == None) and (opts.curve_filename == None):
        print "MODELING MODE WITH NO SHAPE DESCRIPTOR!!!"
    
#---------------get config file name-------------------------------------
    if opts.config != None:   
        inp.config = opts.config
    elif opts.config == None:
        inp.config = ""
        #raise InputError("No config file name was given, please provide one!")
        
#---------------get output name-------------------------------------
    if opts.output != None:   
        inp.outname = opts.output
    elif opts.output == None:
        #if no output folder is provided create pyryresults folder in actual localization
        inp.outname = "pyryresults"
        
#------------write trajectory file -----------------------------------
    if opts.traflfile != None:   
        inp.traflfile = opts.traflfile
    elif opts.traflfile == None:
        print "Data will not be added to trajectory file!!"
        
#------------fullatom file -----------------------------------
    if opts.fullatom != None:   
        inp.fullatom_file = opts.fullatom
    elif opts.fullatom == None:
        inp.fullatom_file = ""
        
#------------save history of moves to output file-------------------------
    if opts.movehistory != None:   
        inp.movehistory_file = opts.movehistory
    elif opts.movehistory == None:
        inp.movehist_file = ""
        
#------------print penalties plot -------------------------
    if opts.score_plot != None:   
        inp.scoreplot = opts.score_plot
    elif opts.score_plot == None:
        inp.scoreplot = ""
        
#-------------wrong option name ------------------------------------
    else: raise InputError("inproper option, try pyry3d.py --help")


#----------- allow auto module ----------------------------
    if opts.auto != None:
        
        if opts.map_filename == None and opts.saxs_filename == None and opts.curve_filename == None and opts.restr_filename == None:
            raise InputError("to run simulations you must provide either complex shape or restraints")
        
        inp.autofolder = opts.auto
        if opts.config == None: InConfig().generate_pyry_inconfig("config_file.txt")    
        instr = InStructures()
        instr.generate_pyry_instructures(str(opts.auto), str(opts.auto)+"_out")
        
        structs = Structures()
        path = str(opts.auto)+"_out.tar"
        inp.structures = structs.get_structures(path, "ig")
        modelsets.append(PDBDirfile(str(opts.auto)+"_out.tar", structs.pdb_files))
        modelset=CompositeModelSet(modelsets)
        
        InSequences().generate_pyry_insequences(str(opts.auto)+".fasta", instr.structures)
        seqs = Sequences()
        inp.sequences = seqs.get_seqs_data(str(opts.auto)+".fasta")
        
        
    else: pass    
        #here check whether restraints or map are provided --> if not do not run simulations or run them with warning that this is pointless.

    if opts.seq_filename == None and opts.auto == None:   
        InSequences().generate_pyry_insequences("sequences.fasta", inp.structures, mode = "ig")
        seqs = Sequences()
        inp.sequences = seqs.get_seqs_data("sequences.fasta")
    return opts



