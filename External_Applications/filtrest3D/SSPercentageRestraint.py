# see license.txt

""" 
Module comment 
SSPercentageRestraint -- class checking percentage ss restraints
"""
import PModel  #getting access to stride_dict
from Restraint import Restraint

class SSPercentageRestraint(Restraint):
    """
    This class stores information about lower and upper percentage bound
    of existing strides of specified type.    
    """
    type=["Secondary structure","Percentage"]
    MAX_DEVIATION=20 # maximum deviation in SS percentage points
    def __init__(self, strideSSType, lowerBound, upperBound, weight=1.0, id=""):
        """
	strideSSType - type of stride to be checked ('H' for helix, 'E' for strand)
	lowerBound - lower percentage bound of strides of selected type
	upperBound - upper percentage bound of strides of selected type
	weight - optional weight
	"""
        self.lowerBound=lowerBound 
        self.upperBound=upperBound
	self.stype=strideSSType
        self.weight=weight
        self.id=id
    
    def autoname(self):
       """ Label for a column with results of the restraint """
       return str(self.lowerBound) + "<=SA<=" + str(self.upperBound)
   
    def check(self, pmodel):
        """Check if SS content satisfies the restraint.
        """
        fraction = self.value(pmodel)
        #now  counting penalty    
        if fraction  <= self.upperBound and fraction >= self.lowerBound: 
            return 0 #number of strides is between given bounds - no penalty
        else:
            deviation = min(abs(self.lowerBound-fraction), abs(fraction-self.upperBound))
        if deviation>self.MAX_DEVIATION:
          deviation = self.MAX_DEVIATION
        return deviation/self.MAX_DEVIATION*self.PENALTY_MAX*self.weight

    def value(self, pmodel):
        """
        This method counts the number of strides of the selected type (H or B),
        and returns penalty  according to given lower and upper percentage bound (times optional weight).
        """
        import sys
        try: 
            test = pmodel.model.get_list()[0].get_list()[0].secondary_structure
        except AttributeError:
            pmodel.stride()
        found = 0.0 #stores the number of strides of specified type in pmodel.stride_dict
        all   = 0
        for chain in pmodel.get_iterator(): #counting number of strides of specified type
            for residue in chain.get_iterator():
                if hasattr(residue, "secondary_structure"):
                    all += 1
                    if residue.secondary_structure in self.stype:
                        found += 1.0
        fraction = (found*100.0/all)
        print >>sys.stderr, "Number of %s is %s" % (self.stype, fraction)
        return fraction
        
