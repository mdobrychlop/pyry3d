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

from random import choice, sample, uniform
from Modules.Simul.Scaler import *
from copy import deepcopy
from Modules.Error.Errors          import InputError

class Mutator(object):
    def __init__(self, scaler):
        self.scaler = scaler
        self.saveHistory = False
    
    def setSaveHistory(self, saveHistory):
        self.saveHistory = saveHistory
        
    def chooseRandomComponent(self, complex):
        componentIndex = choice(complex.movable)
        component = complex.components[componentIndex]
        #component_copy = deepcopy(component)
        #complex.replace_component(componentIndex, component_copy)
        return componentIndex, component
    #
    #def CopyAllComponents(self, complex):
    #    
    #    for i in range(0, len(complex.components)):
    #        component_copy = deepcopy(complex.components[i])
    #        complex.replace_component(i, component_copy)
        
    def mutate(self, complex):
        raise NotImplementedError("Needs subclass implementation of mutate()")
        
    def __str__(self):
        return "Mutator " + self.__class__.__name__
        
class Factory():
    def __init__(self, configuration, simulationState, valueScaler = None, rangeScaler = None, transScaler = None):
        self.configuration = configuration
        self.simulationState = simulationState
        self.scaler = valueScaler if valueScaler != None else self.__prepareDefaultScaler()
        self.rangeScaler = rangeScaler if rangeScaler != None else self.__prepareDefaultRangeScaler()
        self.translationScaler = transScaler if transScaler != None else self.__prepareTranslationRangeScaler()
        self.producer_map = {
            "translation"       : "_Factory__produceTranslation", 
            "translation_all"   : "_Factory__produceTranslationAll", 
            "rotation"          : "_Factory__produceRotation", 
            "rotation_cov"      : "_Factory__produceRotationCov", 
            "rotation_all"      : "_Factory__produceRotationAll", 
            "rotation_whole"    : "_Factory__produceRotationWhole", 
            "exchange"          : "_Factory__produceExchange",
            "exchangeandsample" : "_Factory__produceExchangeAndSample",
            "SimulateDisorder"  : "_Factory__produceSimulateDisorder", 
        }
        self.frequencyCutoffs = {}
        self.frequencySum = 0.
        for frequencyKey in self.configuration.mutation_frequencies.keys():
            if self.configuration.mutation_frequencies[frequencyKey] > 0.:
                self.frequencySum += self.configuration.mutation_frequencies[frequencyKey]
                self.frequencyCutoffs.update({self.frequencySum: frequencyKey})
        
    def produce(self):
        draw = uniform(0., self.frequencySum)
        mutation = ""
        for cutoff in sorted( self.frequencyCutoffs.keys() ):
            if draw < cutoff:
                mutation = self.frequencyCutoffs[cutoff]
                break
        produce = getattr(self, self.producer_map[mutation])
        mutator = produce()
        mutator.setSaveHistory(True)
        return mutator
    
    def __prepareDefaultScaler(self):
        if self.configuration.param_scaling == "on":
            lWhereAreWe = float(self.simulationState.mCurrentIteration) / float(self.simulationState.mIterationCount)
            lRangeIndex = 0
            for lIndex in range(0, len(self.configuration.param_scaling_ranges)):
                if lWhereAreWe < self.configuration.param_scaling_ranges[lIndex]:
                    break
                lRangeIndex = lIndex
            lScales = self.configuration.scaling_ranges[lRangeIndex]
    #        print "Scale ranges:", lScales
            return RandomSignScaler( RandomFloatScaler(None, lScales[0] / 100., lScales[1] / 100.) )
        else:
            return RandomFloatScaler(None, 0., 1.)
        
    def __prepareDefaultRangeScaler(self):
        if self.configuration.param_scaling == "on":
            lWhereAreWe = float(self.simulationState.mCurrentIteration) / float(self.simulationState.mIterationCount)
            lRangeIndex = 0
            for lIndex in range(0, len(self.configuration.param_scaling_ranges)):
                if lWhereAreWe < self.configuration.param_scaling_ranges[lIndex]:
                    break
                lRangeIndex = lIndex
            lScales = self.configuration.scaling_ranges[lRangeIndex]
            #print "Scale ranges:", lScales
            return RangeScaler(lScales[0] / 100., lScales[1] / 100.)
        else:
            return RangeScaler(0., 1.)
        
    def __prepareTranslationRangeScaler(self):
        if self.configuration.param_scaling == "on":
            lWhereAreWe = float(self.simulationState.mCurrentIteration) / float(self.simulationState.mIterationCount)
            lRangeIndex = 0
            for lIndex in range(0, len(self.configuration.param_scaling_ranges)):
                if lWhereAreWe < self.configuration.param_scaling_ranges[lIndex]:
                    break
                lRangeIndex = lIndex
            lScales = self.configuration.scaling_ranges[lRangeIndex]
            #print "Scale ranges:", lScales
            return RangeScaler(0., lScales[1] / 100.)
        else:
            return RangeScaler(0., 1.)
        
    def __produceExchange(self):
        print "Producing exchange"
        mutator = Exchange(self.scaler)
        return mutator
    
    def __produceExchangeAndSample(self):
        print "Producing exchange and sample"
        mutator = ExchangeAndSample(self.rangeScaler, self.configuration, self.simulationState)
        return mutator

    def __produceTranslation(self):
        print "Producing translate"
        mutator = Translation(self.translationScaler, self.configuration.max_trans_vec)
        return mutator

    def __produceTranslationAll(self):
        print "Producing translate all"
        mutator = TranslationAll(self.translationScaler, self.configuration.max_trans_vec)
        return mutator

    def __produceRotation(self):
        print "Producing rotate"
        mutator = Rotation(self.rangeScaler, self.configuration.max_rot_angle)
        return mutator

    def __produceRotationAll(self):
        print "Producing rotate all"
        mutator = RotationAll(self.rangeScaler, self.configuration.max_rot_angle)
        return mutator

    def __produceRotationWhole(self):
        print "Producing rotate whole"
        mutator = RotationWhole(self.rangeScaler, self.configuration.max_rot_angle)
        return mutator

    def __produceRotationCov(self):
        print "Producing rotate covalent"
        mutator = RotationCovalent(self.rangeScaler, self.configuration.max_rot_angle)
        return mutator

    def __produceSimulateDisorder(self):
        print "Producing simulate disorder"
        mutator = SimulateDisorder(self.scaler)
        return mutator

class Translation(Mutator):
    def __init__(self, scaler, maxTranslationVector):
        super(Translation, self).__init__(scaler)
        self.maxTranslationVector = maxTranslationVector
        
    def mutate(self, complex):
        #get random component
        componentIndex, component = self.chooseRandomComponent(complex)
        self.doTranslation(complex, componentIndex, component)
        
    def doTranslation(self, complex, componentIndex, component):
        lVector = []
        
        for i in range(0, len(component.moves.trans_ranges)):
            lRange = [
                max(component.moves.trans_ranges[i][0], component.moves.trans_ranges_sum[i][0]),
                min(component.moves.trans_ranges[i][1], component.moves.trans_ranges_sum[i][1])
            ]
            lVector.append(self.scaler.scale(lRange))
        
        #print "Got vector "+str(lVector)
        
        component = self.performTranslation(complex, componentIndex, lVector, self.saveHistory)
        
        if component.covalent_bonds:
            already_mutated = [componentIndex]
            for covbond in component.covalent_bonds:
                for comp_index in covbond.chains_indexes:
                    if comp_index not in already_mutated:
                        self.performTranslation(complex, comp_index, lVector, self.saveHistory)
                        already_mutated.append(comp_index)        

        return lVector
    
    def performTranslation(self, complex, comp_index, lVector, history):
        comp = complex.components[comp_index]
        component_copy = deepcopy(comp)
        complex.replace_component(comp_index, component_copy)  
        component_copy.translate(lVector, history)
                    #
        if comp_index not in complex.changes:
            complex.add_simul_change_by_index(comp_index)
        if component_copy.moves.limited == True:
            component_copy.check_translation_limits()
        return component_copy

class TranslationAll(Translation):
    def mutate(self, complex):
        #self.CopyAllComponents(complex)
        lVector = []
        lowerBounds = [[], [], []]
        upperBounds = [[], [], []]
        for component_index in complex.movable: #iterate over components
            component = complex.components[component_index]
            if component.moves:
                for i in range(0, len(component.moves.trans_ranges)): # iterate over axi
                    lowerBounds[i].append(max(component.moves.trans_ranges[i][0], component.moves.trans_ranges_sum[i][0]))
                    upperBounds[i].append(min(component.moves.trans_ranges[i][1], component.moves.trans_ranges_sum[i][1]))
                
        
        for axis in range(0, len(lowerBounds)):
            #print "Axis ", axis, " lowerBounds: ", lowerBounds[axis]
            axisRange = [
                max(lowerBounds[axis]),
                min(upperBounds[axis])
                        ]
            #print "Axis ", axis, " range: ", axisRange
            lVector.append(self.scaler.scale(axisRange))
            
        for componentIndex in complex.movable:
            self.performTranslation(complex, componentIndex, lVector, self.saveHistory)
    

class Rotation(Mutator):
    def __init__(self, scaler, maxRotationAngle):
        super(Rotation, self).__init__(scaler)
        self.maxRotationAngle = maxRotationAngle
        self.axis = None
        self.angle = None
    
    def setAxis(self, axis):
        self.axis = axis
    
    def setAngle(self, angle):
        self.angle = angle
    
    def mutate(self, complex):
        componentIndex, component = self.chooseRandomComponent(complex)
        self.doRotation(complex, componentIndex, component)
    
    def setupAxisAndAngle(self, component):
        if self.axis == None:
            lAxi = component.moves.rot_axis
            #for components can not be rotated but their translations are allowed
            try: self.axis = choice(lAxi)
            except: self.axis = "X"
            
        #choose rotation angle
        if self.angle == None:
            if component.moves:
				toScale = [
					max(component.moves.rot_ranges[self.axis][0], component.moves.rot_ranges_sum[self.axis][0]),
					min(component.moves.rot_ranges[self.axis][1], component.moves.rot_ranges_sum[self.axis][1])
				]
            else:
                toScale = self.maxRotationAngle
            
            self.angle = self.scaler.scale(toScale)
            #print "Scaling from", toScale, "default?", self.maxRotationAngle, "chosen angle", self.angle
            
    def doRotation(self, complex, componentIndex, component):
        self.setupAxisAndAngle(component)
        component = self.performRotation(complex, componentIndex)
        
        if component.covalent_bonds:
            rotation_center = -component.mass_centre #zmiana
            
            already_mutated = [componentIndex]
            for covbond in component.covalent_bonds:
                for comp_index in covbond.chains_indexes:
                    if comp_index not in already_mutated:
                        self.performRotation(complex, comp_index, rotation_center)  #zmiana
                        already_mutated.append(comp_index)        
    
        return self.angle, self.axis, [], 0

            
    def performRotation(self, complex, comp_index, rotpoint=[]):
        comp = complex.components[comp_index]

        component_copy = deepcopy(comp)
        complex.replace_component(comp_index, component_copy)
        if len(rotpoint) != 0:
            component_copy.rotate("rotate_whole", self.angle, self.axis, self.saveHistory, rotpoint)
        else:
            component_copy.rotate("", self.angle, self.axis, self.saveHistory)
                
        if comp_index not in complex.changes:
            complex.add_simul_change_by_index(comp_index)  
        #if a component has limited move ranges, rot_ranges must be rescaled
        if component_copy.moves.limited == True:
            component_copy.check_rotation_limits(self.axis)
	return component_copy

class RotationCovalent(Rotation):
    def __init__(self, scaler, maxRotationAngle):
        super(RotationCovalent, self).__init__(scaler, maxRotationAngle)
        self.inpoints = []
        self.cov_index = None
    
    def setInpoints(self, inpoints):
        self.inpoints = inpoints
    
    def setCovIndex(self, covIndex):
        self.cov_index = covIndex
        
    def setupAxisAndAngle(self, component):
        self.axis = "L"
            
        super(RotationCovalent, self).setupAxisAndAngle(component)

    def doRotation(self, complex, componentIndex, component):
        self.setupAxisAndAngle(component)
        component_copy = deepcopy(component)
#by JMK:
        complex.replace_component(componentIndex, component_copy)
        self.inpoints, self.cov_index = component_copy.rotate_around_covalent_bond(self.angle, self.saveHistory)
        if self.cov_index == None: return 1
        for comp_index in component.covalent_bonds[self.cov_index].chains_indexes:
            comp = complex.components[comp_index]
            component_copy = deepcopy(comp)
            complex.replace_component(comp_index, component_copy)
            component_copy.rotate_around_covalent_bond(self.angle, self.saveHistory, self.inpoints, self.cov_index)
                
    
class RotationAll(Rotation):
    def mutate(self, complex):
        #self.CopyAllComponents(complex)
        self.setupAxisAndAngle(complex)
        for componentIndex in complex.movable:
            self.doRotation(complex, componentIndex, complex.components[componentIndex])
#dodane..    
    def doRotation(self, complex, componentIndex, component):
        self.setupAxisAndAngle(component)
        self.performRotation(complex, componentIndex)
    
    def setupAxisAndAngle(self, complex):
        if self.axis == None:
            self.axis = choice(["X", "Y", "Z"])
        
        if self.angle == None:
            self.angle = self.scaler.scale(self.chooseAngleRange(complex, self.axis))

    def chooseAngleRange(self, complex, axis):
        lowerBounds = []
        upperBounds = []
        for component_index in complex.movable:
            component = complex.components[component_index]
            if component.moves:
                lowerBounds.append(max(component.moves.rot_ranges[axis][0], component.moves.rot_ranges_sum[axis][0]))
                upperBounds.append(max(component.moves.rot_ranges[axis][1], component.moves.rot_ranges_sum[axis][1]))
        chosenRange = [max(lowerBounds), min(upperBounds)]
        return chosenRange

class RotationWhole(RotationAll):

#not private any more?
    def performRotation(self, complex, comp_index):
        comp = complex.components[comp_index]
        component_copy = deepcopy(comp)
        complex.replace_component(comp_index, component_copy)
        component_copy.rotate("rotate_whole", self.angle, self.axis, self.saveHistory, complex.centre_of_mass)         
        
        if comp_index not in complex.changes:
            complex.add_simul_change_by_index(comp_index)  
        #if a component has limited move ranges, rot_ranges must be rescaled
        if component_copy.moves.limited == True:
            component_copy.check_rotation_limits(self.axis)
            
class Exchange(Mutator):
    def mutate(self, complex):
        if (len(complex.free) >= 2):
            chosen = sample(complex.free, 2)
            chosen1_copy = deepcopy(complex.components[chosen[0]])
            chosen2_copy = deepcopy(complex.components[chosen[1]])
            complex.replace_component(chosen[0], chosen1_copy)
            complex.replace_component(chosen[1], chosen2_copy)
            
            complex.exchange_components(chosen1_copy, chosen2_copy, self.saveHistory)
            complex.add_simul_change_by_index(chosen[0])  
            complex.add_simul_change_by_index(chosen[1])
            
class ExchangeAndSample(Mutator):
    def __init__(self, scaler, configuration, simulation):
        super(ExchangeAndSample, self).__init__(scaler)
        self.scaler = scaler
        self.maxTranslationVector = configuration.max_trans_vec #maxTranslationVector
        self.maxRotationAngle     = configuration.max_rot_angle #maxRotationAngle
        self.simulation           = simulation
        
        
    def mutate(self, complex ): #, mIterationCount, mCurrentIteration):
        mIterationCount, mCurrentIteration = self.simulation.mIterationCount, self.simulation.mCurrentIteration
        if (len(complex.free) >= 2):
            chosen = sample(complex.free, 2)
            chosen1_copy = deepcopy(complex.components[chosen[0]])
            chosen2_copy = deepcopy(complex.components[chosen[1]])
            complex.replace_component(chosen[0], chosen1_copy)
            complex.replace_component(chosen[1], chosen2_copy)
            
            complex.exchange_components(chosen1_copy, chosen2_copy, self.saveHistory)
            complex.add_simul_change_by_index(chosen[0])  
            complex.add_simul_change_by_index(chosen[1])
            self.sample(complex, mIterationCount, mCurrentIteration, chosen)
        
    def sample(self, complex, mIterationCount, mCurrentIteration, chosen_indexes):
        
	self.simulation.calculateRating(complex)
	best_score = complex.simulation_score
	
	mut_nr = (mIterationCount - mCurrentIteration)*0.01 #from simul class
	if mut_nr > 100:
	    mut_nr = 100
	elif mut_nr < 10: mut_nr = 10 
	
	#print "*****", mut_nr, self.mCurrentIteration, self.mIterationCount
	for i in range(0, int(mut_nr)):
	    mCurrentIteration += 1 #from simul class
	    #self.mIterationCount -= 1
	    #print "!!!", self.mCurrentIteration, self.mIterationCount
	    print "searching..", i
	    complex.changes = []
	    NewComplex = complex.deepcopy_score() #deepcopy_score?
	    

	    componentIndexNew = sample([chosen_indexes[0], chosen_indexes[1]], 1)[0]
	    mutationNew = sample(["Rotation", "Translation"], 1)[0]
            #MutatedComponentNew = NewComplex.components[componentIndexNew]
	    MutatedComponentNew = deepcopy(NewComplex.components[componentIndexNew])
	    #self.calculateRating(NewComplex)
	    #self.calculateRating(aComplex)
	    if mutationNew == "Rotation":
                mutator = Rotation(self.scaler, self.maxRotationAngle)
                mutator.doRotation(NewComplex, componentIndexNew, MutatedComponentNew)
	    else:
                mutator = Translation(self.scaler, self.maxTranslationVector)
                mutator.doTranslation(NewComplex, componentIndexNew, MutatedComponentNew)
	    self.simulation.calculateRating(NewComplex)

	    if NewComplex.simulation_score > best_score:
                complex.add_simul_change_by_index(componentIndexNew)
                complex.replace_component(componentIndexNew, MutatedComponentNew)
                self.simulation.calculateRating(complex)
                best_score = NewComplex.simulation_score
		
	

	#print "check whether exchange was proper", self.calculateRating(aComplex)

class SimulateDisorder(Mutator):
    def mutate(self, complex):
        if len(complex.with_disorders) == 0: raise InputError("You cannot apply Simulate Disorder mutation for components with no flexible/disordered fragments")
        componentIndex = choice(complex.with_disorders)
        component = complex.components[componentIndex]
	
	component_copy = deepcopy(component)
        complex.replace_component(componentIndex, component_copy)
	
        #simulate disorder
        ddstruct = component.simulate_disorder(component_copy.pyrystruct.struct, "simul")
        if ddstruct:
            if componentIndex not in complex.changes:
                complex.add_simul_change_by_index(componentIndex)
#            self.mutated_components[componentIndex] = component 
