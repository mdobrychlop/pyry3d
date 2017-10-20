#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Modules.Trans.Interaction import *

class RestraintFilter (object):
    def __init__(self, changes = []):
        self.changes = changes
    
    def __call__(self, interactions):
        return self.filter(interactions)

    def filter(self, interactions):
        result_restraints = []
        if len(self.changes):
            for restraint in interactions:
                for change in self.changes:
                    if change in restraint.chains:
                        result_restraints.append(restraint)
                        break
        else:
            result_restraints = interactions
        return result_restraints
        
class RegularRestraintFilter ( RestraintFilter ):
    def filter(self, interactions):
        interactions = super(RegularRestraintFilter, self).filter(interactions)
        result_restraints = []
        for restraint in interactions:
            if not isinstance(restraint, SymmetryDistances):
                result_restraints.append(restraint)
        return result_restraints
        
class SymmetryRestraintFilter ( RestraintFilter ):
    def filter(self, interactions):
        interactions = super(SymmetryRestraintFilter, self).filter(interactions)
        result_restraints = []
        for restraint in interactions:
            if isinstance(restraint, SymmetryDistances):
                result_restraints.append(restraint)
        return result_restraints
