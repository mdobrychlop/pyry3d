#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# 
# www.genesilico.pl 
#

__author__ = "Joanna & Dominik Kasprzak"
__copyright__ = "Copyright 2010, The PyRy3D Project"
__credits__ = ["Janusz Bujnicki"]
__license__ = "GPL"
__version__ = "0.1.0"
__maintainer__ = "Joanna & Dominik Kasprzak"
__email__ = "jkasp@amu.edu.pl"
__status__ = "Prototype"

from random import uniform, random

class Scaler(object):
    def __init__(self, other = None):
        self.other = other
    
    def preScale(self, value):
        return value if self.other == None else self.other.scale(value)
        
    def scale(self, value):
        raise NotImplementedError("Needs subclass implementation of scale()")


class InvariantScaler(Scaler):
    def scale(self, value):
        return self.preScale(value)
        
class RandomFloatScaler(Scaler):
    def __init__(self, other = None, minimumPercentage = 0.0, maximumPercentage = 100.0):
        super(RandomFloatScaler, self).__init__(other)
        self.minimumPercentage = minimumPercentage
        self.maximumPercentage = maximumPercentage
        
    def scale(self, value):
        value = self.preScale(value)
        return uniform(value * self.minimumPercentage, value * self.maximumPercentage)

class RandomSignScaler(Scaler):
    def scale(self, value):
        value = self.preScale(value)
        return value if uniform(0, 1) > 0.5 else value * -1.0
        
class AddConstantScaler(Scaler):
    def __init__(self, constant, other = None):
        super(AddConstantScaler, self).__init__(other)
        self.constant = constant
        
    def scale(self, value):
        value = self.preScale(value)
        return value + self.constant

class AddPercentScaler(Scaler):
    def __init__(self, percent, other = None):
        super(AddPercentScaler, self).__init__(other)
        self.percent = percent
        
    def scale(self, value):
        value = self.preScale(value)
        return value + self.percent * value 
       
class MultiplyScaler(Scaler):
    def __init__(self, constant, other = None):
        super(MultiplyScaler, self).__init__(other)
        self.constant = constant
        
    def scale(self, value):
        value = self.preScale(value)
        return value * self.constant
    
class RangeScaler(Scaler):
    def __init__(self, minimumPercentage = 0.0, maximumPercentage = 100.0, other = None):
        super(RangeScaler, self).__init__(other)
        self.minimumPercentage = minimumPercentage
        self.maximumPercentage = maximumPercentage
        
    def scale(self, value):
        aMinimumValue = value[0]
        aMinimumValueSign = aMinimumValue / abs(aMinimumValue) if aMinimumValue else 1
        aMaximumValue = value[1]
        aMaximumValueSign = aMaximumValue / abs(aMaximumValue) if aMaximumValue else 1
        isSameSign = (aMinimumValue * aMaximumValue) >= 0
        divisionPoint = (aMinimumValue + aMaximumValue) / 2 if isSameSign else 0
        lowerRange = [
            ((aMinimumValue - divisionPoint*aMinimumValueSign) * self.minimumPercentage) + divisionPoint*aMinimumValueSign,
            ((aMinimumValue - divisionPoint*aMinimumValueSign) * self.maximumPercentage) + divisionPoint*aMinimumValueSign
        ]
        upperRange = [
            ((aMaximumValue - divisionPoint*aMaximumValueSign) * self.minimumPercentage) + divisionPoint*aMaximumValueSign,
            ((aMaximumValue - divisionPoint*aMaximumValueSign) * self.maximumPercentage) + divisionPoint*aMaximumValueSign
        ]
        lScaled = uniform(lowerRange[0], lowerRange[1]) if random() < 0.5 else uniform(upperRange[0], upperRange[1])
        #print "Scaled from "+str(aMinimumValue)+" to "+str(aMaximumValue)+" to "+str(lScaled)
        return lScaled


