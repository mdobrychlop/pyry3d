from Modules.Trans.RestraintChecker import *
from Modules.Trans.RestraintScorer import *
from Modules.Trans.ScoreSaver import *
from Modules.Trans.Interaction import *

class RestraintScorerFactory(object):
    def __init__(self, complex):
        self.complex = complex
        
    def createScorer(self, restraint, savers = None):
        if isinstance(restraint, SymmetryDistances):
            return self.__getSymmetryRestraintScorer(restraint, savers)
        if isinstance(restraint, DoubleComponentInteraction):
            return self.__getPairRestraintScorer(restraint, savers)
        if isinstance(restraint, SingleComponentInteraction):
            return self.__getSingleRestraintScorer(restraint, savers)
        if isinstance(restraint, LogicalRestraint):
            return self.__getLogicalRestraintScorer(restraint, savers)

    def __getSymmetryRestraintScorer(self, restraint, savers = None):
        restraints_components = []
        componentA = self.complex.get_component_by_chain(restraint.chains[0])
        componentB = self.complex.get_component_by_chain(restraint.chains[1])
        restraints_components.append(componentA.pyrystruct)
        restraints_components.append(componentB.pyrystruct)
        #calc = SymmetricRestraintChecker(componentA.pyrystruct, componentB.pyrystruct)
        calc = SymmetricRestraintChecker(restraints_components)
        scorer = SingleRestraintScorer(restraint, calc, self.complex, savers)
        return scorer

    def __getPairRestraintScorer(self, restraint, savers = None):
        restraints_components = []
        componentA = self.complex.get_component_by_chain(restraint.chains[0])
        componentB = self.complex.get_component_by_chain(restraint.chains[1])
        restraints_components.append(componentA.pyrystruct)
        restraints_components.append(componentB.pyrystruct)
        if restraint.type == "relation":
            componentC = self.complex.get_component_by_chain(restraint.chains[2])
            componentD = self.complex.get_component_by_chain(restraint.chains[3])
            restraints_components.append(componentC.pyrystruct)
            restraints_components.append(componentD.pyrystruct)

        calc = PairRestraintChecker(restraints_components)
        scorer = SingleRestraintScorer(restraint, calc, self.complex, savers)
        return scorer

    def __getSingleRestraintScorer(self, restraint, savers = None):
        component = self.complex.get_component_by_chain(restraint.chains[0])
        calc = SingleComponentRestraintChecker(component)
        scorer = SingleRestraintScorer(restraint, calc, self.complex, savers)
        return scorer
        
    def __getLogicalRestraintScorer(self, restraint, savers = None):
        scorer = LogicalRestraintScorer(restraint, self, savers)
        return scorer
        

