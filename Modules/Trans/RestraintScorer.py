#class RestraintScorer(object):
#    def __init__(self, filter, checker, scoreSaver):
#        self.filter = filter
#        self.checker = checker
#        self.scoreSaver = scoreSaver
#    
#    def score(self, restraints):
#        targetRestraints = self.filter(restraints)
#        score = self.checker.check(targetRestraints)
#        self.scoreSaver( score )

class SingleRestraintScorer(object):
    def __init__(self, restraint, scoreCalculator, complex, savers = None):
        self.restraint = restraint
        self.calculator = scoreCalculator
        self.complex = complex
        self.savers = savers
    
    def score(self):
        score = self.calculator.calculate(self.restraint)
        self.restraint.set_score(score)
        if (self.savers):
            for saver in self.savers:
                saver(self.restraint.uid, score)
                
class LogicalRestraintScorer(object):
    def __init__(self, restraint, factory, savers = None):
        self.restraint = restraint
        self.factory = factory
        self.savers = savers
    
    def score(self):
        self.restraint.update_score(self.factory)
        score = self.restraint.get_score()
        if (self.savers):
            for saver in self.savers:
                saver(self.restraint.uid, score)
