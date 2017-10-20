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

"""
module stores all errors classes used in a program
"""

class PyryError(Exception): pass

class InputError(PyryError): pass
class SequenceError(InputError) : pass
class DensityMapError(InputError): pass
class StructureError(InputError) : pass
class FiltrestError(InputError): pass
class RestraintError(InputError): pass
class PyRyStructureError(InputError): pass

class ConfigError(PyryError): pass

class TransError(PyryError): pass
class ComponentError(TransError): pass
class Component_movesError(TransError) : pass
class InteractionError(TransError): pass
class RestraintCheckerError(TransError): pass
class Complex_mapError(TransError): pass
class ComponentRepresentationError(TransError): pass
class ShapeError(TransError) : pass


class SimulError(PyryError): pass
class PyRyComplexError(SimulError): pass


class ClusterError(PyryError): pass

class FilterError(PyryError): pass

class TrajectError(PyryError): pass

class Out_setError(PyryError): pass

class Final_modelError(PyryError): pass

class EvaluateError(PyryError): pass


    
    
    
    
    
    
    
    
    
    
    
    
    
