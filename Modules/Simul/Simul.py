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

from math                           import ceil, log, pow, fabs, exp
from random                         import randint, uniform, sample, random
from copy                           import deepcopy, copy
from sys							import maxint, exit
from numpy import percentile

from Modules.Error.Errors           import SimulError #to represent this module errors
#usage: raise SimulError("No structure folder given, I have nothing to do!")

from Modules.Simul.Scaler import *
from Modules.Simul.Mutator          import *

from Modules.Trans.Component        import Component
#class representing one component = structure 

from Modules.Simul.Complex          import PyRyComplex
from Modules.Simul.PlotComplex      import PlotComplex
from Modules.Simul.Simul_cpp		import cpp_anneal_complex, cpp_calculate_rating, \
										cpp_choose_result_save_procedure, \
										cpp_metropolis_compare, cpp_save_last_simulated_complex, \
										cpp_init, cpp_save_simulation_results, transfer_to_pyrycomplex
from Modules.Simul.Reductions_cpp	import *
from Modules.Constans.Logfile       import logfile, movehist
from Modules.Simul.cpp				import pyry3d_cpp
from Bio.PDB.Structure              import Structure
from Bio.PDB.Model                  import Model
from Bio.PDB                        import PDBIO

from Modules.Constans.Logfile       import logfile


class PlotComplex(object):
    
    def __init__(self, complex, iteration):
	self.sim_score = complex.simulation_score 
	self.clashes   = complex.clashes     * complex.clashes_penalty
	self.restraints = complex.restraints * complex.restraints_penalty
	self.outbox    = complex.outbox      * complex.outbox_penalty
	self.mapfill   = complex.freespace   * complex.freespace_penalty
	self.densities = complex.density     * complex.density_penalty
	self.step_nr   = iteration

class Simul(object):
    #Mutations = ["Rotate", "Translate", "Exchange", "SimulateDisorder"]
    class SimulationMethods:
        Genetic = 1
        SimulatedAnnealing = 2
        ReplicaExchange = 3
    
    class ReductionMethods:
        Roulette = 1
        Tournament = 2
        Cutoff = 3
        
    def __init__(self, configuration):
        self.configuration = configuration
        #   default values for simulation parameters
        self.mSimulationMethod = Simul.SimulationMethods.SimulatedAnnealing
        self.mGeneticReductionMethod = Simul.ReductionMethods.Roulette
        
        self.mOutputCount = 1
        self.mIterationCount = 5
        self.mCurrentIteration = 0
        self.mMaxPoolSize = 5
        self.mMaxTemperature = 1000
        self.mPScalRangeBounds = []
        self.mPScalValues = []
        self.mReplicaExchangeSteps = 1
        self.mReplicaTemperatures = []
        self.mReplicaOdd = 0
        
        self.rejected = 0   #to calculate number of rejected complexes!
     
        self.mInteractions = None
        self.mDensityMap = None



#--------- simulation score weights -------------------------------
        self.collision_penalty  = 1.
        self.restraints_penalty = 1.
        self.freespace_penalty  = 1.
        self.outbox_penalty     = 1.
        self.density_penalty    = 1.
        self.symmetry_penalty   = 1.    #should be used for optimalization only (for symmetrical complexes)

#################
        self.x =0  #check number of accepted changes; used only for testing purposes 
####################
        
        self.mRejectionThreshold = 0.6  ##  Currently not used, see reducePopulation() for details
        self.mMaximumPopulationCount = 0 ## Setting this to a number greater than
        #0 will case population to be reduced to the mMaximumPopulationCount best specimens each evolution run.
        
        self.mPool = list() ## Initialize a pool for the simulation with an empty list
        
        self.scores = []   #here unredundant simul scores are written
            
        self.traflfile         = None #trajectory file if provided
	self.saved_iterations  = []   #to store indexes of complexes already saved on disc
	self.complexes_to_save = []   #if user defines steps of simulation to save a list stores these complexes here
	self.saved_complexes   = []   #list of complexes scores
	self.plot              = False #does the user want to get score plot?
	self.mutated_components = {}
	self.lBestSimulation_score =  float(-(2**1024 - 2**971))
	# minimum possible value of float
	#-1000000000000000000000000000000000000000000000000000 :D
	self.lBestScoreComplex    = None # best score from whole simulation
	self.lBestScoreComplex_niter = None #best in last N steps
	
	self.Reheating            = False
	self.RejectionsFreq       = None
	
	self.server				  = False # if true, program runs simulation in c++

	
    #   Simulated annealing methods
    def __annealComplex (self, aComplex, aIteration, aMaxTemperature):
	
        lNewComplex = aComplex.deepcopy_score()
        lNewComplex.set_iteration_index(aIteration)
        self.performMutation(lNewComplex)
        self.calculateRating(lNewComplex)
	lNewComplex.calculate_centre_of_complex()
	lCurrentTemp = self.__temperature(aIteration, aMaxTemperature)
        lNewComplex.temp = lCurrentTemp
        
	#self.save_simulation_results(aIteration, [lNewComplex], "actual")
	    
	if ( self.__metropolisCompare(lNewComplex, aComplex, lCurrentTemp) ):#Accept the change
            del aComplex
	    #print "!!Accepted!!!"
	    
	    #print "mut components", self.mutated_components.keys()
	#    if self.simul_params.movehistory:    
	#	for comp_index in xrange(0,len(lNewComplex.components)):
	#	    if comp_index not in self.mutated_components.keys():
	#		
	#		lNewComplex.components[comp_index].rot_history.append([[0,0,0],[0]])
	#		lNewComplex.components[comp_index].trans_history.append([0,0,0])
			#print "not changed?", comp_index, lNewComplex.components[comp_index].trans_history, len(lNewComplex.components[comp_index].trans_history)
	    
	    
            return lNewComplex
        else:
            del lNewComplex
            #print "!!Rejected!!!"
            return aComplex
    
    def anneal (self):
	"""
	Performs simulated annealing procedure; Returns True id change was
	accepted and False when it was discard

	If anything will be changed here, you have to check cpp_anneal_complex
	in Simul_cpp.py file, and check if c++ version still works correctly
	"""
        lOldScore = self.mPool[0].simulation_score if not self.server \
				else pyry3d_cpp.get_simulation_score(0)
	# regular annealing

	if (self.mLastAnnealingRestart == None or self.mLastAnnealingRestart == 0 or self.Reheating == False):
		if self.server:
			cpp_anneal_complex(0, pyry3d_cpp.get_step(), self.mMaxTemperature, self)
		else:
			self.mPool[0] = self.__annealComplex(self.mPool[0], self.mCurrentIteration, self.mMaxTemperature)
	# reheated annealing
	else:
		if self.server:
			cpp_anneal_complex(0, pyry3d_cpp.get_step()-self.mLastAnnealingRestart, self.mMaxTemperature/2, self)
		else:
			self.mPool[0] = self.__annealComplex(self.mPool[0], \
					self.mCurrentIteration-self.mLastAnnealingRestart, self.mMaxTemperature/2)
			
	sim_score = self.mPool[0].simulation_score if not self.server else \
		pyry3d_cpp.get_simulation_score(0)
	
	if (lOldScore != sim_score): #True if score changed
	    self.mAnnealingRejections = 0
	    return True
	else:
	    self.mAnnealingRejections += 1
	    if (self.mAnnealingRejections > self.RejectionsFreq): #TODO Needs parametrization
		self.mLastAnnealingRestart = self.mCurrentIteration if not self.server \
			else pyry3d_cpp.get_step()
		self.mAnnealingRejections = 0
	    return False
	
    def calculateRating (self, aComplex):
        """
	Calculates new simulation score for a specified structure (aComplex).
	"""
	aComplex.calculate_simulation_score(self.mInteractions, self.mDensityMap, self.simul_params)	
	
    def __calculateScaledWeight(self, scoring_element):
        """
        calculates scaled value for a given simulation score element and returns
        new value
        """
        range = scoring_element[0] - scoring_element[1]
        new_penalty = self.scaleParameter(fabs(range), True)
        if range < 0:
            new_penalty = scoring_element[1] - new_penalty
        else:
            new_penalty = scoring_element[1] + new_penalty
        return new_penalty
    
    def collect_complexes_to_save(self, i, aComplex):
	"""
	collects complexes which will be returned to the user according
	to defined simulation steps to save and number of best outcomplexes
	"""
	
	if self.simul_params.out_steps[0] <= i <= self.simul_params.out_steps[1]:
	    if len(self.complexes_to_save) < self.simul_params.struct_nr:
		self.complexes_to_save.append(deepcopy(aComplex))
	    else:
		self.complexes_to_save.sort( key=lambda Complex: Complex.simulation_score)
		for el in self.complexes_to_save:
		    if aComplex.simulation_score > el.simulation_score:
			self.complexes_to_save.remove(el)
			self.complexes_to_save.append(deepcopy(aComplex))
			break
	
    def __cutoffReduction (self):
	"""
	reduced population size based on cutoff value
	"""
        lNewPool = self.mPool
        if (len(self.mPool) <= self.mMaxPoolSize):
            return lNewPool
        lCutoffIndex = self.mMaxPoolSize
        if lCutoffIndex > 0:
            lSelectedIndices = range(0, lCutoffIndex)
            lNewPool = self.mPool[:lCutoffIndex]
        lRejectedSpecimens = []
        for lIndex in range(0, len(self.mPool)-1):
            if lSelectedIndices.count(lIndex) == 0:
                lRejectedSpecimens.append(self.mPool[lIndex])
        for lRejectedComplex in lRejectedSpecimens:
            del lRejectedComplex
        return lNewPool
    
    #   Methods related to genetic simulation
    def geneticSimulation (self):
	"""
	Performs genetic simulation procedure

	Any changes here migth change behaviour of c++ version
	"""
        # Let's iterate over a copy of the gene pool (we're modifying the pool, so making a copy is necessary)
        print "Starting iteration"
        
        mPoolSize = len(self.mPool) if not self.server else pyry3d_cpp.pool_size()	
        for lComplex in xrange(0, mPoolSize):
            if self.server:
                pyry3d_cpp.set_pool_size(pyry3d_cpp.pool_size() + 1)
                new_index = pyry3d_cpp.perform_mutation(lComplex)
                cpp_calculate_rating(new_index, self)
            else:
                lNewComplex = deepcopy(self.mPool[lComplex])
                self.mPool.append(lNewComplex)     
                self.performMutation(lNewComplex)
                self.calculateRating(lNewComplex)
		lNewComplex.calculate_centre_of_complex()
         
        print "New generation created"
        # A new generation was created, let's reduce the population. Only the strong or lucky may survive!
        print "Sorting"
        if self.server:
			pyry3d_cpp.sort_by_simulation_score()
        else:
			self.mPool.sort(key=lambda Complex: Complex.simulation_score, reverse=True)
        print "Reducing"
        # w srodku self.server!!!
        self.__reducePopulation()
        print "Re-sorting"
        if self.server:
			pyry3d_cpp.sort_by_simulation_score()
        else:
			self.mPool.sort(key=lambda Complex: Complex.simulation_score, reverse=True)
        print "Survived:"
        if self.server:
            for lComplex in xrange(0, pyry3d_cpp.pool_size()):
				print str(pyry3d_cpp.get_simulation_score(lComplex))
        else:
            for lComplex in self.mPool:
				print lComplex.simulation_score
	    
    def get_plot_values(self):
	return self.saved_complexes
	
    def getResult (self):
        """Returns simulation results.
        If result count is smaller than specimen pool size, self.mOutputCount
	best structures are returned. Otherwise all specimens are given."""
        if (len(self.mPool) > self.mOutputCount):
            lResult = self.mPool[0:self.mOutputCount]
        else:
            lResult = self.mPool
        return lResult
    
    def performMutation (self, aComplex):
        """Performs a mutation on a specified structure (aComplex).
        Mutation type is random and can be a rotation or translation.
        Mutation type parameters are random as well (see @see setMaxTranslationDistance
	and @see setMaxRotationAngle for details).
        @todo Add support for another type of mutation - a component is replaced
	by a component from another structure."""
	
        #for first complex all components should be treated as changed!!
        aComplex.clean_simul_changes()
	self.__scaleWeights(aComplex)
        factory = Factory(self.configuration, self)
        mutator = factory.produce()
        mutator.mutate(aComplex)
        return
		
    def __metropolisCompare (self, aComplexA, aComplexB, aCurrentTemp = 0):
        if (aCurrentTemp == 0):
            aCurrentTemp = self.mMaxTemperature
        
        lDiff = aComplexA.simulation_score - aComplexB.simulation_score
	#asia:
	if aCurrentTemp == 0. : aCurrentTemp = 0.0001
	##
        if ( (-lDiff)/aCurrentTemp < 700.0): #There is a limit to math.exp at approx 700.0, after that it throws OverflowError
            lAcceptanceProbability = 1.0 / (1.0 + exp( (-lDiff)/aCurrentTemp ) )
        else:
            lAcceptanceProbability = 0.0
        lRandom = random()
        
        #print "Diff:", lDiff, "temp:", aCurrentTemp, "acc prob:", lAcceptanceProbability, "random:", lRandom , "scores:", aComplexA.simulation_score, aComplexB.simulation_score
        
        if ( (lDiff > 0) or (lRandom < lAcceptanceProbability) ):    #Accept the change
            print "CHANGE ACCEPTED!!!"
            self.x += 1
	    
            return True
        return False
    
    def populatePool (self):
        """Private method used at the start of simulation.
        Populates specimen pool with copies of the original structure, filling
	25% of it (rounded up)
	"""
        lInitialPoolSize = ceil(0.25 * self.mMaxPoolSize)
        print "Populating with initial", lInitialPoolSize, "clones"
        lComplex = self.mPool[0]
        while (len(self.mPool) < lInitialPoolSize):
            lClonedComplex = deepcopy (lComplex)
            self.mPool.append (lClonedComplex)
    
    def __reducePopulation (self):
        """
	Performed after each iteration, reduces the population to fit
	self.mMaxPoolSize specimens.

	Note that we call c++ functions here
        """
        if not self.server:
            if (len(self.mPool) <= self.mMaxPoolSize):
                return
        else:
        	if (pyry3d_cpp.pool_size() <= self.mMaxPoolSize):
        		return
        if (self.mGeneticReductionMethod == Simul.ReductionMethods.Cutoff):
            if self.server:
                cpp_cutoff(self)
            else:
                self.mPool = self.__cutoffReduction()
        elif (self.mGeneticReductionMethod == Simul.ReductionMethods.Roulette):
            if self.server:
                cpp_roulette(self)
            else:
                self.mPool = self.__rouletteReduction()
        elif (self.mGeneticReductionMethod == Simul.ReductionMethods.Tournament):
            if self.server:
                cpp_tournament(self)
            else:
                self.mPool = self.__tournamentReduction()
        else:
            raise SimulError ("Wrong population reduction method selected")
	    
    def replica_exchange(self):
        """
	Performs replica exchange procedure

	If anything changed here, you have to actualize c++ version
         """
        if self.server:
            if (pyry3d_cpp.pool_size() != len(self.mReplicaTemperatures)):
                raise "Replica exchange: Pool size should match temperature count"
        else:
            if (len(self.mPool) != len(self.mReplicaTemperatures)):
                raise "Replica exchange: Pool size should match temperature count"
        
        for x in range(0, len(self.mReplicaTemperatures)):
            if self.server:
                pyry3d_cpp.set_reptemp(x, self.mReplicaTemperatures[x])
                cpp_anneal_complex(x, self.mCurrentIteration, self.mReplicaTemperatures[x], self);
            else:
				self.mPool[x].reptemp = self.mReplicaTemperatures[x]
				self.mPool[x] = self.__annealComplex(self.mPool[x], self.mCurrentIteration, self.mReplicaTemperatures[x]);
	  
        if ((self.mCurrentIteration % self.mReplicaExchangeSteps) == 0):
			
            if self.server:
				for i in range(self.mReplicaOdd, pyry3d_cpp.pool_size() - 1, 2):
					if(cpp_metropolis_compare(self.mMaxTemperature, i, i+1, \
							self, self.mReplicaTemperatures[i])):
						print "Exchanging replicas", i, "and", i+1
						pyry3d_cpp.exchange(i, i+1)
            else:
				for i in range(self.mReplicaOdd, len(self.mPool) - 1, 2):
					if (self.__metropolisCompare(self.mPool[i], self.mPool[i+1], \
							self.mReplicaTemperatures[i])):
						print "Exchanging replicas", i, "and", i+1
						lTempComplex = self.mPool[i]
						self.mPool[i] = self.mPool[i+1]
						self.mPool[i+1] = lTempComplex
						#exchange
            self.mReplicaOdd = (self.mReplicaOdd + 1) % 2
	
    def __rouletteReduction (self):
	"""
	performs reduction of population size based on roulette criteria
	"""
        lMinScore = abs( self.mPool[0].simulation_score )
        lMaxScore = abs( self.mPool[len(self.mPool)-1].simulation_score )
        lDiff = lMaxScore - lMinScore
        lTemperature = self.mMaxTemperature * log(self.mIterationCount - self.mCurrentIteration + 1, self.mIterationCount)

        med_array = []

        lScoreSum = 0
        for lComplex in self.mPool:
            med_array.append(lComplex.simulation_score)

        quantile = percentile( med_array, 25)

        for lComplex in self.mPool:
        	#calculates sum of deviations from median for complexes better than
        	#middle one
        	if lComplex.simulation_score > quantile:
        		lScoreSum += lComplex.simulation_score - quantile
  
        lRouletteWheel = dict()
        lCurrentProb = 0.0
        for lComplexIndex in xrange(0, len(self.mPool)): #Calculate probability for each complex, add it to "The wheel"
            lComplex = self.mPool[lComplexIndex]
            lProb = 0 if (lComplex.simulation_score < quantile or lScoreSum == 0) else \
            	-(quantile - lComplex.simulation_score) / (lScoreSum)

            lCurrentProb += lProb
            lRouletteWheel[ lComplexIndex ] = lCurrentProb
            #print "Score:", abs(lComplex.simulation_score), ", prob:", lProb
        #The probabilities will sum up to 1.
            
        #print "Wheel of fortune stakes:"
        #print lRouletteWheel
        lNewPool = []
        lSelectedIndices = []
        while len(lNewPool) < self.mMaxPoolSize: #Draw from the wheel until new pool is full
            lSelectedIndex = -1
            lDrawResult = random() * lCurrentProb
            #print "Draw result:", lDrawResult
            for i in xrange(0, len(lRouletteWheel)):
                if lRouletteWheel[i] > lDrawResult:
                    lSelectedIndex = i
                    break
            if (lSelectedIndex > -1) and (lSelectedIndices.count(lSelectedIndex) == 0):
                lNewPool.append(self.mPool[lSelectedIndex])
                lSelectedIndices.append(lSelectedIndex)
            #    print "Selected index:", lSelectedIndex
            #break
        lRejectedSpecimens = []
        for lIndex in range(0, len(self.mPool)-1):
            if lSelectedIndices.count(lIndex) == 0:
                lRejectedSpecimens.append(self.mPool[lIndex])
        for lRejectedComplex in lRejectedSpecimens:
            del lRejectedComplex
        return lNewPool
    
    def save_simulation_results(self, i, lBestScoreComplex, fname=False):
	"""
	Saves results of simulation to trafl file (if requested by the user) and
	to pdb files
	
	results are saved every N steps (defined by the user) and best scored
	complex in last N iterations is saved
	
	start/stop number of iterations to save (requested by the user)
	"""	
	if self.traflfile and not self.server: # and (i not in self.saved_iterations):
	    for cx in lBestScoreComplex: cx.add_trajectory(i, cx.temp)
		
	for replica in lBestScoreComplex: 
	    if replica == None:
	    	continue  
	    if replica.simulation_score not in self.scores:
		self.scores.append(replica.simulation_score)
		if self.mReplicaTemperatures: replica.save_pdb(i,replica.reptemp, fname) #self.mReplicaTemperatures[j])
		else: replica.save_pdb(i, replica.temp, fname)
		if not self.server:
			logfile.write_file(\
			    "cx score for iteration "+str(i+1)+" is\t"+str(round(replica.simulation_score,3))\
			    +"\tcomponents: restraints: "+str(round(replica.restraints,3) * replica.restraints_penalty)+" "+\
			    "collisions: "+str(round(replica.clashes,3)    * replica.clashes_penalty)+" "+\
			    "map filling: "+str(round(replica.freespace,3) * replica.freespace_penalty)+" "+\
			    "outbox atoms: "+str(round(replica.outbox,3)     * replica.outbox_penalty)+" "+\
			    "density filling: "+str(round(replica.density,3)  * replica.density_penalty)+" "+\
			    "\tACTUAL WEIGHTS are respectively: "+str(replica.restraints_penalty)+", "+\
			    str(replica.clashes_penalty)+", "+str(replica.freespace_penalty)+", "+\
			    str(replica.outbox_penalty)+", "+str(replica.density_penalty)+"\n")

		if self.plot == True:
		    plot_cx = PlotComplex(replica, i)
		    self.saved_complexes.append(plot_cx)
	    
	    else: logfile.write_file("Rejected\t"+str(replica.simulation_score)+"\titer_nr\t"+str(i+1)+"\n")
		
##@TODO: put into separate function
	#if self.simul_params.movehistory:
	#    for replica in lBestScoreComplex: 
	#	movehist.write_message("COMPLEX \t"+str(i)+ "\t"+str(replica.simulation_score)+"\t"+str(replica.temp))
	#	for component in replica.components:
	#	    movehist.write_message(str(component.pyrystruct.chain)+"\tROTATIONS\t"+str(component.rot_history)\
	#			+"\tTRANSLATIONS\t"+str(list(component.trans_history)))
		    
	self.saved_iterations.append(i)
		
    def __scaleComponentsNumber(self, aComplex):
	"""
	Chooses number of components to be mutated during simulation
	When scaling is on than scaling procedure is used to calculate this number
	otherwise random number is selected
	"""
        if self.scaling == "on":
            lMutatedComponentsNumber = int(round(self.scaleParameter(len(aComplex.movable)-1, True), 0)) + 1
        elif self.scaling == "off":
	    if len(aComplex.movable) == 1:
		lMutatedComponentsNumber = 1
	    else:
	        lMutatedComponentsNumber = randint(1, len(aComplex.movable))
        return lMutatedComponentsNumber
    
    def scaleParameter (self, aMaximumValue, aOnlyPositive = False):
	"""
	scales parameter based on user defined scaling ranges and returns its
	new value
	"""
        lWhereAreWe = float(self.mCurrentIteration) / float(self.mIterationCount)
        lRange = 0
        for lIndex in range(0, len(self.mPScalRangeBounds)):
            #print "Range", lIndex, "value:", self.mPScalRangeBounds[lIndex]
            if lWhereAreWe < self.mPScalRangeBounds[lIndex]:
                break
            lRange = lIndex
        #print "Where are we? We're at", lWhereAreWe, ", range index", lRange
        lScales = self.mPScalValues[lRange]
        #print "Scale ranges:", lScales
        lMultiplier = uniform(lScales[0], lScales[1])
        lScaledParam = aMaximumValue * lMultiplier
        if (not aOnlyPositive) and (uniform(0,1) > 0.5):
            lScaledParam = lScaledParam * -1.0
        #print "Scaled value:", lScaledParam
        return lScaledParam

    def scaleParameterRange (self, aMinimumValue, aMaximumValue):
	"""
	scales parameter from range aMinimumValue - aMaximumValue based on user
	defined scaling ranges and returns its  new value
	"""
        lWhereAreWe = float(self.mCurrentIteration) / float(self.mIterationCount)
        lRange = 0
        for lIndex in range(0, len(self.mPScalRangeBounds)):
            #print "Range", lIndex, "value:", self.mPScalRangeBounds[lIndex]
            if lWhereAreWe < self.mPScalRangeBounds[lIndex]:
                break
            lRange = lIndex
        #print "Where are we? We're at", lWhereAreWe, ", range index", lRange
        lScales = self.mPScalValues[lRange]
        #print "Scale ranges:", lScales
        lMultiplier = uniform(lScales[0], lScales[1])
        if (aMinimumValue > aMaximumValue):
            lTemp = aMinimumValue
            aMinimumValue = aMaximumValue
            aMaximumValue = lTemp
        lMiddle = (aMinimumValue + aMaximumValue) / 2.0
        lMax = aMaximumValue - lMiddle
        lScaledParam = lMax * lMultiplier
        if (uniform(0,1) > 0.5):
            lScaledParam = lScaledParam * -1.0
        lScaledParam += lMiddle
        ####
        #print "Scaled value:", lScaledParam, " Middle: ",lMiddle
        return lScaledParam
        
    def __scaleWeights(self, aComplex):
        """
           for each simulation score element weight is scaled according to
           simulation progress and ranges defined by the user
        """
        if self.scaling == "on":
            new_out_box_penalty    = self.__calculateScaledWeight(self.simul_params.outbox_penalty)
            new_free_space_penalty = self.__calculateScaledWeight(self.simul_params.freespace_penalty)
            new_collision_penalty  = self.__calculateScaledWeight(self.simul_params.clashes_penalty)
            new_restraints_penalty = self.__calculateScaledWeight(self.simul_params.restraints_penalty)
            new_density_penalty    = self.__calculateScaledWeight(self.simul_params.density_penalty)
            new_symmetry_penalty   = self.__calculateScaledWeight(self.simul_params.symmetry_penalty) 
            new_chi2_penalty       = self.__calculateScaledWeight(self.simul_params.chi2_penalty)
            new_rg_penalty         = self.__calculateScaledWeight(self.simul_params.rg_penalty)
	    
            aComplex.set_penalties_weights(new_collision_penalty, new_restraints_penalty, \
					new_free_space_penalty, new_out_box_penalty, \
					new_density_penalty, new_symmetry_penalty, new_chi2_penalty, new_rg_penalty)
        else:
            aComplex.set_penalties_weights(self.simul_params.clashes_penalty[0],
                       self.simul_params.restraints_penalty[0],
                       self.simul_params.freespace_penalty[0],
                       self.simul_params.outbox_penalty[0],
                       self.simul_params.density_penalty[0],
                       self.simul_params.symmetry_penalty[0],
                       self.simul_params.chi2_penalty[0],
                       self.simul_params.rg_penalty[0])
    
    def setDensityMap (self, aMap):
        """Sets density map used for score calculation"""
        self.mDensityMap = aMap
	
    def setIndicators(self, simul_params):
        self.simul_params = simul_params
	
    def setIterationCount (self, aCount):
        """Sets number of performed simulation iterations"""

        self.mIterationCount = aCount
	
    def setInteractions (self, aInteractions):
        """Sets interaction objects used for score calculation"""
        self.mInteractions = aInteractions
	
    def setMaximumPoolSize (self, aPoolSize):
        """Sets maximum size of population pool.
        A larger pool size allows more specimen variety but requires more memory
        and computation power to process. It is advised to experiment with the best settings.
        The pool size cannot be smaller then output structure count.
        """
        if (aPoolSize >= self.mOutputCount):
            self.mMaxPoolSize = aPoolSize
        else:
            self.mMaxPoolSize = self.mOutputCount
    
    def setMutationFrequencies(self, rotations, rotations_cov, translations, exchange, exchangeandsample, simul_dd, translations_all, \
    		rotations_all, rotations_whole):
        """ sets frequencies of all mutations' types"""
        self.rotation_freq = rotations
	self.rotation_cov_freq = rotations_cov
        self.translation_freq = translations
        self.rotation_whole_freq = rotations_whole
        self.exchange_freq = exchange
	self.exchangeandsample = exchangeandsample
        self.simul_dd = simul_dd
        self.translation_all_freq = translations_all
	self.rotation_all_freq = rotations_all
	
        #print "SIMUL mut params", self.rotation_freq, self.translation_freq, \
        #self.exchange_freq, self.simul_dd, self.rotation_cov_freq, self.rotation_all_freq, self.translation_all_freq

    def setParameterScalingBoundries (self, aBoundaries):
        """ aBoundries - an array of scaling group upper boundries."""
        for lRange in aBoundaries:
            if lRange > 1:
                lRange = lRange / 100.0
            self.mPScalRangeBounds.append(lRange)
        print aBoundaries
        
    def setParameterScalingValues (self, aValues):
        for lRange in aValues:
            for i in range(0, len(lRange)):
                if lRange[i] > 1:
                    lRange[i] = lRange[i] / 100.0
            self.mPScalValues.append(lRange)
        #print aValues
	
    def setReheatingParams(self, reheat, rejections):
        self.Reheating = reheat
	self.RejectionsFreq = rejections
	
    def setScoresPlot(self):
	self.plot = True
	
    def setReductionMethod (self, aMethod):
        if aMethod == "Roulette": aMethod = 1
        elif aMethod == "Tournament": aMethod =2
        elif aMethod == "Cutoff": aMethod =3
        else: raise SimulError ("Wrong reduction method selected")
        self.mGeneticReductionMethod = aMethod
	
    def setRejectionThreshold (self, aThreshold):
        """Score rejection threshold. Not used ATM"""
        self.mRejectionThreshold = aThreshold
	
    def setReplicaExchangeSteps(self, exchange_freq):
	self.mReplicaExchangeSteps = exchange_freq
	
    def setReplicaTemperatures(self, rep_temp):
	self.mReplicaTemperatures = rep_temp
	
    def setResultCount (self, aCount):
        """Sets number of structures returned at the end of simulation"""
        self.mOutputCount = aCount
        if (self.mMaxPoolSize < aCount):
            self.mMaxPoolSize = aCount
	    
    def setScalingOption(self, option):
        "sets scaling to on or off"
        self.scaling = option
	 
    def setServer(self, optimized):
		self.server = optimized
	    
    def setSimulationMethod (self, aMethod):
        if aMethod == "simulatedannealing": aMethod = 2
        elif aMethod == "genetic": aMethod =1
        elif aMethod == "replicaexchange": aMethod = 3
        else: raise SimulError ("Wrong simulation method selected")
        self.mSimulationMethod = aMethod
	
    def setStartingComplex (self, aComplex):
        """Sets the initial structure."""
        self.mPool = []
        self.mPool.append (aComplex)
	
    def setStepsToSave(self, traj_step):
        """Sets which steps will be added to trajectory"""
        self.mTrajStep = traj_step
        
    def setTemperature (self, aTemperature):
        self.mMaxTemperature = aTemperature
        
    def setTraflfile(self, traflfile):
        """ Sets tajectory file if given by the user """
        self.traflfile = traflfile

    def start(self):
        """Starts the simulation process.

        If anything has been changed in python's data structure you have to:
        	1. Change this in c++ version
        	2. Initialize this data structure in c++
        	3. Perhaps deallocate this data structure in c++ (if there were allocations)

        Also, be very careful if changes here happen

        @todo is it true that only one component is changed in single simul step??"""
        print "Simulation started..\nComplex Score\n"
        print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t"%("Total", "Restraints", \
				"Collisions", "MapFill", "Outbox", "DensityFill", "Symmetry")
	
        if self.server:
            self.mPool[0].interactions = self.mInteractions
	    cpp_init(self)

        if (self.mSimulationMethod == Simul.SimulationMethods.Genetic):
            if self.server:
                pyry3d_cpp.populate_pool(self.mMaxPoolSize)
            else:
				self.populatePool ()
        if (self.mSimulationMethod == Simul.SimulationMethods.SimulatedAnnealing):
            self.mLastAnnealingRestart = 0
            self.mAnnealingRejections = 0
        if (self.mSimulationMethod == Simul.SimulationMethods.ReplicaExchange):
            for i in range(1, len(self.mReplicaTemperatures)):
                if self.server:
                    pyry3d_cpp.deepcopy_complex(-1)
                else:
                    self.mPool.append(self.mPool[0])
            if self.server:
                pyry3d_cpp.set_pool_size(len(self.mReplicaTemperatures))
        
        if self.simul_params.movehistory:
            if self.server:
                pyry3d_cpp.deepcopy_complex_to_best_score_complex(0)
            else:
                self.lBestScoreComplex = deepcopy(self.mPool[0]) #deepcopy()
        else: 
            self.lBestScoreComplex = None
        if self.simul_params.save_res_mode == "niter":
            if self.server:
                pyry3d_cpp.deepcopy_complex_to_best_score_niter_complex(0)
        else:
                self.lBestScoreComplex_niter = deepcopy(self.mPool[0]) #deepcopy()
	
		#set best score counter    
        self.lBestSimulation_score = self.mPool[0].simulation_score

        for i in xrange(0, self.mIterationCount) :
			self.mCurrentIteration = i
			self.simul_params.iter_nr = i
			if self.server:
				print "Iteration " + str(i) + ", pool size:" + \
						str(pyry3d_cpp.pool_size())
			else:
				print "Iteration " + str(i) + ", pool size:" + str(len(self.mPool))
		   
			if (self.mSimulationMethod == Simul.SimulationMethods.Genetic):
				self.geneticSimulation()
			elif (self.mSimulationMethod == Simul.SimulationMethods.SimulatedAnnealing):
				self.anneal()
			elif (self.mSimulationMethod == Simul.SimulationMethods.ReplicaExchange):
				self.replica_exchange()
			else:
				raise SimulError ("Wrong simulation method selected")
			self.mutated_components = {}
			
			if self.server:
				cpp_choose_result_save_procedure(i, self)
			else:
				self.__choose_result_save_procedure(i)

		#save last complex!!
        if self.server:
			cpp_save_last_simulated_complex(i, self, Simul.SimulationMethods)
        else:
			self.__save_last_simulated_complex(i)
        print "Simulation params:", self.mOutputCount, self.mIterationCount, \
				self.mMaxTemperature, self.mPScalRangeBounds, self.mPScalValues, \
				self.mSimulationMethod, self.mGeneticReductionMethod
			

        print "number of accepted changes during simulation", self.x

        #deallocation
        pyry3d_cpp.free_cpp()

        if self.simul_params.movehistory and not self.server:
			print "best complex is..", self.lBestScoreComplex.iteration_index, self.lBestScoreComplex.simulation_score
			return self.lBestScoreComplex, self.mPool[0]
        return 0, 0

	
    def __save_last_simulated_complex(self, i):
	"""
	i - iteration index
	"""
	if self.simul_params.save_res_mode == "outsteps":
	    self.complexes_to_save.sort( key=lambda Complex: Complex.iteration_index)
	    for outcomplex in self.complexes_to_save:
	        self.save_simulation_results(outcomplex.iteration_index, [outcomplex])
		print "should be saved", outcomplex.iteration_index
		
	if self.simul_params.save_res_mode == "niter":
	    if (self.mSimulationMethod == Simul.SimulationMethods.ReplicaExchange): self.save_simulation_results(i, self.mPool)
	    else: self.save_simulation_results(i, [self.mPool[0]])        
        
	
    def __choose_result_save_procedure(self, iteration_number):
	"""
	depending on type of results saving defined by the user from the following options:
	eachbetter - save on disc each complex with score better than for previous strucutres
	niter      - if user wants to save results after each N steps of simulation
	outsteps   - when a user defines numbers of steps which will be saved on disc.
	
	also the method copies the complex if its score is the best from already obtained
	in order to store coordinates of best complex and at the end of run to recreate a fullatom model from reduced representation.
	"""
	    
	#save each new complex with better score
	if self.simul_params.save_res_mode == "eachbetter":
	    if self.mPool[0].simulation_score > self.lBestSimulation_score:
		self.save_simulation_results(iteration_number, [self.mPool[0]])
		self.lBestSimulation_score = self.mPool[0].simulation_score
	    #elif (iteration_number % self.mTrajStep == 0):
	#	self.save_simulation_results(iteration_number, [self.mPool[0]])
		
	#save complex with best score to retrieve full atom model
	if self.simul_params.movehistory:
	    if (self.mPool[0].simulation_score > self.lBestScoreComplex.simulation_score):
		#print "move histoey saving best complex", self.mPool[0].simulation_score, self.lBestScoreComplex.simulation_score
		self.lBestScoreComplex = deepcopy(self.mPool[0]) 
		    
	#save best complex from last N iterations
	if self.simul_params.save_res_mode == "niter":
	    if (iteration_number % self.mTrajStep == 0):
		if (self.mSimulationMethod == Simul.SimulationMethods.ReplicaExchange): self.save_simulation_results(iteration_number, self.mPool)
		self.save_simulation_results(iteration_number, [self.lBestScoreComplex_niter])
		#print "zapis!, kompleksy//", self.lBestScoreComplex_niter.simulation_score, self.mPool[0].simulation_score
		self.lBestScoreComplex_niter = deepcopy(self.mPool[0])
	    else:
		if self.mPool[0].simulation_score > self.lBestScoreComplex_niter.simulation_score:
		    #print "ZAMIANA!!!", self.mPool[0].simulation_score, self.lBestScoreComplex_niter.simulation_score
		    self.lBestScoreComplex_niter = deepcopy(self.mPool[0])
		    
	#save only if iteration if user want complex generated in this simulation step to be saved	    
	if self.simul_params.save_res_mode == "outsteps":
	    self.collect_complexes_to_save(iteration_number, self.mPool[0])
	    
	#print "actual complexes rozmanzany komplex", self.mPool[0].simulation_score, self.mPool[0].iteration_index
	#print "actual bests:", lBestScoreComplex.simulation_score, lBestScoreComplex_niter.simulation_score
	
    def __temperature (self, step, aTemp = 0):
	"""
	returns new temperature value
	"""
        if (aTemp == 0):
            aTemp = self.mMaxTemperature
        alpha = 0.999
        return aTemp * pow(alpha, step)
	
    def tournamentReduction (self):
	"""
	reduces population size based on tournament criteria
	"""
        lNewPool = []
        lSelectedIndices = []
        
        while (len(lNewPool) < self.mMaxPoolSize):
            lContender1 = -1
            lContender2 = -1
            while ( (lContender1 == -1) or (lContender2 == -1) ):
                if (lContender1 == -1):
                    lContender1 = randint(0, len(self.mPool) - 1)
                    if (lSelectedIndices.count(lContender1) > 0):
                        lContender1 = -1
                if (lContender2 == -1):
                    lContender2 = randint(0, len(self.mPool) - 1)
                    if (lSelectedIndices.count(lContender2) > 0) or (lContender2 == lContender1):
                        lContender2 = -1
            #print "Fight between", lContender1, self.mPool[lContender1].simulation_score, "and", lContender2, self.mPool[lContender2].simulation_score
            if (self.mPool[lContender1].simulation_score > self.mPool[lContender2].simulation_score):
                lNewPool.append(self.mPool[lContender1])
                lSelectedIndices.append(lContender1)
            else:
                lNewPool.append(self.mPool[lContender2])
                lSelectedIndices.append(lContender2)
        #print "Winners:", lSelectedIndices
        lRejectedSpecimens = []
        for lIndex in range(0, len(self.mPool)-1):
            if lSelectedIndices.count(lIndex) == 0:
                lRejectedSpecimens.append(self.mPool[lIndex])
        for lRejectedComplex in lRejectedSpecimens:
            del lRejectedComplex
        return lNewPool
    
	
if __name__=='__main__':
    print "Hello!"
