
class PlotComplex(object):
    
    def __init__(self, complex, iteration):
	self.sim_score = complex.simulation_score 
	self.clashes   = complex.clashes     * complex.clashes_penalty
	self.restraints = complex.restraints * complex.restraints_penalty
	self.outbox    = complex.outbox      * complex.outbox_penalty
	self.mapfill   = complex.freespace   * complex.freespace_penalty
	self.densities = complex.density     * complex.density_penalty
	self.step_nr   = iteration
