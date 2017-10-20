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


#math
from math import sqrt


"""
        module containt mathematical operations on vectors
"""

def sqr(x):
        """calculates square """
        return x**2

def count_dist(pt1, pt2):
        """calcuates distance between 2 points in 3D space"""
        #pt1 and pt2 must be LISTS of x,y,z coordinates
        di = sqrt(sqr(pt1[0]-pt2[0])+sqr(pt1[1]-pt2[1])+sqr(pt1[2]-pt2[2]))
        return di

if __name__=='__main__':
    print count_dist([257.561,-52.247, 12.788 ], [191.753,31.53, -40.49 ])
    
