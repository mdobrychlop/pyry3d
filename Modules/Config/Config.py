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


from Modules.Error.Errors           import ConfigError, Component_movesError
from Modules.Constans.Logfile       import logfile
from Modules.Trans.Component        import Movable
import os



#@TODO: clean this module!!!
#@TODO: remove code duplications
#@TODO: check data correctness provided by the user and default values



class Config:
    """
        class which represents simulation configuration data:
        gets simulation parameters and delegates them to SIMUL module;
        gets simulation output data and delegates them to OUT_set module;

        #prepare input data for simulation:
        #outside a map, no restraints kept, collisions, empty spaces in the map
        #give a Simul module information about number os simulation steps
    """
    def __init__(self):
        self.anntemp                  = 10         #annealing temperature
        self.disabled_mutations       = []          #names of mutations with frequencies of 0
        self.clashes_penalty          = [1., 1.]      #penalty for collisions between complex components
        self.curve                    = False       #True if user provides curve file
        self.components_indexes       = {}         #how components are indexed in complex.components list
        self.chi2_penalty             = [0.,0.]     #penalty for disagreement with SAXS curve according to chi2 value
        self.crysol_path              = ""          #path to CRYSOL installed on Ubuntu system
        self.density_penalty          = [0., 0,]    #penalty ranges for occuping points of map with low density
        self.exchange_freq            = 0.0        #frequency of components' exchange mutation
        self.exchangeandsample_freq      = 0.0        #frequency of components' exchange mutation
        self.freespace_penalty       = [1., 1.]    #penalty for empty spaces in density map
        self.identify_disorders       = False        #True if PyRy3D should check for disorders, False if not (faster)
        self.iter_nr                  = 0           #simulation step number
        self.kvol                     = 1.0        #how many complex volumes will define density map
        self.linked_components      = {}          #if some somponents are linked by covalent bond
        self.is_chi2_defined          = False
        self.is_rg_penalty_defined    = False
        self.max_rot_angle            = 5.          #max rotation angle in single simulation move
        self.max_trans_vec            = [5.,5.,5.]  #max translation vector in single simulation move
        self.maxpoolsize              = 100         #maximum pool size for genetic algorithm
        self.movable                  = []          #list of components with defined movements
        self.movehistory              = ""          #filename with history of moves for all saved complexes
        self.niter                    = 0           #which steps will be added to trajectory
        self.outbox_penalty          = [1., 1.]    #penalty for beeing outside simulation box
        self.out_steps                = ["FIRST","LAST"] #which steps to include in output data? default None
        self.param_scaling            = "off"        #scaling procedure can be "on" or "off"
        self.param_scaling_ranges     = [0,25,50]   #values for scaling ranges during simulation
        self.scaling_ranges           = [[50,100], [25,50], [0,25]] #values of scaled params
        self.param_scaling_range1     = [50,100]    #values of scaled params in first range
        self.param_scaling_range2     = [25,50]     #values of scaled params in second range
        self.param_scaling_range3     = [0,25]      #values of scaled params in third range
        self.rg_penalty               = [0.,0.]     #penalty for disagreement with SAXS curve according to Rg value
        self.rg_val                   = 0.0         #Rg value given by user
        self.reductmethod             = "Roulette"  #reduction method for genetic algorithm
        self.reheat                   = False
        self.reheat_rejections        = None
        self.replica_exchange_freq    = 1           #how often replicas will be exchanged; default every 10% of simulation steps
        self.repl_exch_freq           = False       #parameter replicaexhchangefreq not defined by the user
        self.replica_temps            = [400, 350, 300, 250, 200, 150, 100, 50, 25, 0]          #list of replica simulations' temperatures
        self.representation           = "ca"        #structure representation for components; #CA - only calfas/c4'; cg - coarse grain, fa - full atom (default)
        self.restraints_penalty       = [1., 1.]    #penalty for not keeping user defined restraints
        self.mutation_frequencies     = {"rotation" : 0.5, "translation" : 0.5}
        self.rotation_freq            = 0.5        #frequency of component's rotation mutation
        self.rotation_all_freq        = 0.0
        self.rotation_cov_freq        = 0.0
        self.rotation_whole_freq      = 0.0
        self.saxsradius              = 0.0
        self.symmetry_penalty         = [0.,0.]
        self.write_eachbetter         = False
        self.save_res_mode            = "outsteps"  #by default only best energy complex is saved
        self.shapedesc                = True        #True when modeling is with use of map or saxs shape
        self.simbox                   = 1.2         #how many times bigger than map is a simulbox
        self.simboxradius             = 2.          #radius of a single grid cell in Angstroms
        self.simgridtype              = "cubic"     #type of simulation grid; at the moment only cubic grid is available
        self.simmethod                = "simulatedannealing" #simulation algorithm
        self.simul_dd_freq            = 0.0        #frequency of components' simulate disorder mutation
        self.simul_steps              = 100         #number of simulation step to perform
        self.start_orient             = False       #has a user set first complex conformation
        self.struct_nr                = 1           #number of best scores out structures; default 0.1simul_steps
        self.threshold                = None          #min density value accepted in simulation
        self.translation_freq         = 0.5        #frequency of component's translation mutation
        self.translation_all_freq     = 0.0
        self.required_clashes_penalty = False
        self.required_clashes_penalty_allatoms = False
        
        self.is_density_defined         = False
        self.is_freespace_defined       = False
        self.kvol_given                 = False
        self.threshold_given            = False
            
    def __str__(self):
        return "%s %s %s %s %s %s %s" % \
          ( self.simul_steps, self.freespace_penalty, \
           self.clashes_penalty, self.outbox_penalty, self.restraints_penalty,\
           self.out_steps, self.struct_nr)
        
    def __open_config_file(self,filename):
        
        try:
            fh = open(filename, "r")
        except ConfigError: print "Cannot open config file!!", filename
        return fh
    
    def __parse_covalent_links(self, line):
        try: linked_components = eval(line.split()[2].strip()) #.replace("[","").replace("]","").replace(",",""))
        except ConfigError: print "Please provide correct format for COVALENT_BONDS parameter in configuration file"
        #linked_components = line.split()[2].replace("[","").replace("]","").replace(",","")

        
        component = line.split()[1]
        if component not in self.linked_components.keys():
            for el in linked_components:
                if el == component:
                    raise ConfigError("Component %s cannot be bound with %s "%(component, el))
            #self.linked_components[component] = linked_components
            
           
        at1 = eval(line.split()[3])
        at2 = eval(line.split()[4])
        
        covbond = CovalentBond(linked_components, at1, at2)
        if self.linked_components.has_key(component):
            self.linked_components[component].append(covbond)
        else:
            self.linked_components[component] = [covbond]

        #else:
        #    raise ConfigError("You provided more than one covalent bond info for a component:%s"%(component))
        #print "!!!!!", self.linked_components
        
    def __parse_crysol_path(self, line):
        "parse path to crysol"
        path = line.split()[1]
        try:
            os.path.exists(path)
        except ConfigError: print "Please provide valid path to CRYSOL program on your disk"
        self.crysol_path = path
        
    def __parse_reheating_param(self, line):
        
        line = line.split()
        reheat = line[1]
        reheat_rejections = line[2]
        if reheat.upper() not in ["TRUE", "FALSE"]: raise ConfigError("Reheat is of value True or False")
        if reheat_rejections.isalpha(): raise ConfigError("Rejection must be float value")
        elif reheat.upper() == "TRUE" and float(reheat_rejections) > self.simul_steps : raise ConfigError ("Rejection frequency cannot be larger than number of simulation steps")
        elif float(reheat_rejections) <= 0 : raise ConfigError ("Rejection frequency cannot be lower than 0")
        
        if reheat.upper() == "TRUE": self.reheat = True
        else: self.reheat = False
        self.reheat_rejections = float(reheat_rejections)
    
    def __parse_score_weights_penalties(self, line, param_name):
                
        values = [line.split()[1], line.split()[2]]
        if len(values) != 2 : raise ConfigError(str(param_name).upper()+" parameter you must provide 2 values")
        if values[0].upper() == "X" or values[1].upper() == "X": return 1 
        if values[0].isalpha() or values[1].isalpha(): raise ConfigError(str(param_name).upper()+"must be DIGID not alfanumeric value")
        elif float(values[0]) >= 0 and float(values[1]) >=0 : setattr(self, param_name, [float(values[0]), float(values[1])])
        else: raise ConfigError(str(param_name).upper()+"value must be positive value!")
        
        if param_name == "freespace_penalty":
            if float(values[0]) == 0. and float(values[1]) == 0:
                self.is_freespace_defined = False
            else:
                self.is_freespace_defined = True
                
        if param_name == "density_penalty":
            if float(values[0]) == 0. and float(values[1]) == 0:
                self.is_density_defined = False
            else:
                self.is_density_defined = True
                
        if param_name == "chi2_penalty":
            if float(values[0]) == 0. and float(values[1]) == 0:
                self.is_chi2_defined = False
            else:
                self.is_chi2_defined = True
                
        if param_name == "rg_penalty":
            if float(values[0]) == 0. and float(values[1]) == 0:
                self.is_rg_penalty_defined = False
            else:
                self.is_rg_penalty_defined = True
                
        
    def __parse_mutation_frequencies(self, line, param_name):
                
        value = line.split()[1]
        floatValue = float(value)
        if value.upper() == "X": return 1
        #mutation_name = param_name.split("_")[0]
        mutation_name = param_name.replace("_freq","")
        if mutation_name == "simul_dd": mutation_name = "SimulateDisorder"
        if floatValue == 0: self.disabled_mutations.append(mutation_name)
        if value.isalpha() : raise ConfigError(str(param_name)+"must be DIGID not alfanumeric value")
        elif  0. <= floatValue <= 1. : setattr(self, param_name, floatValue)
        else: raise ConfigError("Frequency of"+str(param_name).upper()+"must be defined in range from 0 to 1")
        return {mutation_name: floatValue}

        
    def __parse_positive_params(self, line, param_name, type = "float"):
        """
        """
        value = line.split()[1]
        if value.upper() == "X": return 1
        if value.isalpha() : raise ConfigError(str(param_name)+"must be DIGID not alfanumeric value")
        if not param_name == "threshold":
            try:
                float(value) < 0
            except:
                raise ConfigError(str(param_name).upper()+"must be positive value! Use dots for float numbers.")
        if type == "float":
            setattr(self, param_name, float(value))
        if type == "int":
            setattr(self, param_name, int(value))
                                
    def __parse_param_value(self, line, param_name, list_of_values):
        """
        """
        value = line.split()[1]
        found = False
        
        if value.upper() == "X" : return 1
        if value.isdigit() : raise ConfigError(str(param_name).upper()+"must be alfanumeric not DIGIT value")
        for val in list_of_values:
            if value.lower() == val.lower():
                found = True
                break
        if found == False:
            raise ConfigError(str(param_name).upper()+" has only possible values "+"".join(list_of_values) )
        if param_name == "simmethod" and value.upper() == "X": setattr(self,param_name, "simulatedannealing")
        else: setattr(self,param_name, value.lower())
        
    def __parse_rot_angle(self, line):
        """
        """
        maxrot = line.split()[1]
        if maxrot.upper() == "X": return 1
        if maxrot.isalpha() : raise ConfigError("MAXROT must be DIGID not alfanumeric value")
        if -360.0 > float(maxrot) > 360.: raise ConfigError("Rotation angle cannot be larger than 360 degrees or smaller than -360")
        else: self.max_rot_angle = float(maxrot)
        
    def __parse_outsteps(self, line):
        """
        """
        out_steps = [line.split()[1], line.split()[2]]
        if len(out_steps) != 2 : raise ConfigError("for OUTSTEPS parameter you must provide 2 values")
        if (out_steps[0].upper() == "X" and out_steps[1].upper() == "X") \
        or (out_steps[0].upper() == "FIRST" and out_steps[1].upper() == "LAST"):
            self.out_steps.append(1)
            self.out_steps.append(self.simul_steps-1)
        elif float(out_steps[1]) > self.simul_steps:
            raise ConfigError("Steps value cannot be larger than number of simulations!")
        elif (out_steps[0] > 0 and out_steps[1] > 0): self.out_steps = [int(line.split()[1]), int(line.split()[2])]
        else: raise ConfigError("OutSteps value must be positive value!")
        
    def __parse_trans_vector(self, line):
        """
        """
        trans_vec = [line.split()[1], line.split()[2], line.split()[3]]
        if len(trans_vec) != 3 : raise ConfigError("for MAXTRANS parameter you must provide 3 values")
        if trans_vec[0].upper() == "X" or trans_vec[1].upper() == "X" or trans_vec[2].upper() == "X": return 1
        elif trans_vec[0].isalpha() or trans_vec[1].isalpha() or trans_vec[2].isalpha() : raise ConfigError("MAXTRANS must be DIGIDS not alfanumeric value")
        self.max_trans_vec = [float(trans_vec[0]), float(trans_vec[1]), float(trans_vec[2])]
        
    def __parse_scaling(self, line, param_name, counter):
        """
        """
        if counter == 3:
            values = [float(line.split()[1]), float(line.split()[2]), float(line.split()[3])]
        elif counter == 2:
            values = [float(line.split()[1]), float(line.split()[2])]
        
        percentValues = []
        for val in values:
            if str(val).isalpha(): raise ConfigError(param_name.upper()+"must be DIGITS not alfanumeric value")
            elif 100 < val < 0: raise ConfigError(param_name.upper()+"must be from range 0 to 100")
            percentValues.append(val/100.0)
            
        setattr(self, param_name, percentValues)
        return percentValues
        
    def __parse_movestate(self, line):
        """
        """
        chain_name = line.split()[1]
        move_state = line.split()[2]
        if move_state.lower() == "fixed" or move_state.lower() == "movable": pass
        else: raise ConfigError("unknown state %s"%(move_state))
        move_state = Movable(chain_name, line)
        #print "@@@@@@@@@@2", move_state
        self.movable.append(move_state)
        
    def __parse_kvol(self, line):
        """
        """
        kvol = line.split()[1]
        if kvol.isalpha() : raise ConfigError("KVOL must be DIGIDS not alfanumeric value")
        if kvol.upper() == "X": return 1
        else:
            if (float(kvol) > 10 or float(kvol) < 0):
                raise ConfigError("Volume you provided is not in the range from 0 to 10!%s  "%(kvol))
            elif float(kvol) == 0:
                raise ConfigError("Volume must be larger than 0!")
            self.kvol = float(kvol)
            
    def __parse_bool_param(self,line, param_name):
        val = line.split()[1]
        if val.upper() not in ["TRUE", "FALSE"]:
            raise ConfigError('%s can be "True" or "False"'%(param_name.upper()))
        else:
            if val.upper() == "TRUE": val = True
            else: val=False
            setattr(self, param_name, val)
    
    def __parse_values_list(self, line, param_name):
        """
        retrieves list of digits of unknown length
        """
        li = line.split()[1:]
        temps = []
        for el in li:
            #do not parse elements after comment indicator '#'
            if "#" in el:
                break
            elif el.isdigit(): temps.append(float(el))
        setattr(self, param_name, temps)        
                
    def parse_config_file(self, filename, curvefile, shapedesc = True):
        """
            parses input config file which contains simulation parameters and
            scoring function elements (e.g. penalties for collisions)
        Parameters:
        -----------
            filename   : name of file with simulation parameters
            comp_names : names of complex components
        Returns:
        --------
            config object storing all simulation data
        """
        self.shapedesc = shapedesc
        self.curve = curvefile
        if filename == "":
            self.save_logfile()
            return 1
        
        fh = self.__open_config_file(filename)
        
        #to check whether both params are not defined at the same time
        for line in fh:
            if line.startswith("#"): continue
            
            elif line.upper().startswith("STEPS"):
                self.__parse_positive_params(line,'simul_steps', "int")
                
            elif line.upper().startswith("REHEAT"):
                self.__parse_reheating_param(line)
                
        #parse score weights ---------------------------------    
            elif line.upper().startswith("OUTBOX"):
                self.__parse_score_weights_penalties(line,'outbox_penalty')
                
            elif line.upper().startswith("MAP_FREESPACE"):
                self.__parse_score_weights_penalties(line, 'freespace_penalty')
                
            elif line.upper().startswith("CLASHES "):
                self.required_clashes_penalty = True 
                self.__parse_score_weights_penalties(line, 'clashes_penalty')
                
            elif line.upper().startswith("CLASHES_ALLATOMS"):
                self.required_clashes_penalty_allatoms = True
                self.__parse_score_weights_penalties(line, 'clashes_penalty')
                
            elif line.upper().startswith("RESTRAINTS"):
                self.__parse_score_weights_penalties(line, 'restraints_penalty')
                
            elif line.upper().startswith("DENSITY"):
                self.__parse_score_weights_penalties(line, 'density_penalty')
                
            elif line.upper().startswith("SYMMETRY"):
                self.__parse_score_weights_penalties(line, 'symmetry_penalty')
                
            elif line.upper().startswith("CHI2"):
                self.__parse_score_weights_penalties(line, 'chi2_penalty')
                
            elif line.upper().startswith("RG "):
                self.__parse_score_weights_penalties(line, 'rg_penalty')

        #parse mutation frequencies-----------------------------
            elif line.upper().startswith("ROTATION_FREQ"):
                self.mutation_frequencies.update( self.__parse_mutation_frequencies(line,'rotation_freq') )
                
            elif line.upper().startswith("ROTATION_COV_FREQ"):
                self.mutation_frequencies.update( self.__parse_mutation_frequencies(line,'rotation_cov_freq') )
                
            elif line.upper().startswith("TRANSLATION_FREQ"):
                self.mutation_frequencies.update( self.__parse_mutation_frequencies(line,'translation_freq') )
                
            elif line.upper().startswith("TRANSLATION_ALL_FREQ"):
                self.mutation_frequencies.update( self.__parse_mutation_frequencies(line,'translation_all_freq') )
                
            elif line.upper().startswith("ROTATION_ALL_FREQ"):
                self.mutation_frequencies.update( self.__parse_mutation_frequencies(line,'rotation_all_freq') )
                
            elif line.upper().startswith("ROTATION_WHOLE_FREQ"):
                self.mutation_frequencies.update( self.__parse_mutation_frequencies(line,'rotation_whole_freq') )
                
            elif line.upper().startswith("EXCHANGE_FREQ"):
                self.mutation_frequencies.update( self.__parse_mutation_frequencies(line,'exchange_freq') )
                
            elif line.upper().startswith("EXCHANGESAMPLE_FREQ"):
                self.mutation_frequencies.update( self.__parse_mutation_frequencies(line,'exchangeandsample_freq') )
                
            elif line.upper().startswith("ROTATION_WHOLE_FREQ"):
                self.__parse_mutation_frequencies(line,'rotation_whole_freq')

            elif line.upper().startswith("SIMUL_DD_FREQ"):
                self.mutation_frequencies.update( self.__parse_mutation_frequencies(line,'simul_dd_freq') )
            
            elif line.upper().startswith("WRITE_N_ITER"):
                self.__parse_positive_params(line,'niter', "int")
                self.out_steps = []
            
            elif line.upper().startswith("WRITE_EACHBETTER"):
                self.__parse_bool_param(line, "write_eachbetter")
                self.out_steps = []
                
            elif line.upper().startswith("OUT_STEPS"):
                self.__parse_outsteps(line)
                
            elif line.upper().startswith("STRUCT_NR"):
                self.__parse_positive_params(line,'struct_nr', "int")           
                
            elif line.upper().startswith("MAXROT"):
                self.__parse_rot_angle(line)
                
            elif line.upper().startswith("MAXTRANS"):
                self.__parse_trans_vector(line)
                    
            elif line.upper().startswith("COMPONENT_REPRESENTATION"):          
                self.__parse_param_value(line, "representation", ["fa","ca", "cacb", "3p", "sphere", "ellipsoid"])
                    
            elif line.upper().startswith("SCALEPARAMS"):
                self.__parse_param_value(line, "param_scaling", ["on", "off"])
                
            elif line.upper().startswith("PARAMSCALINGRANGES"):
                self.__parse_scaling(line, "param_scaling_ranges", 3)               
                                
            elif line.upper().startswith("PARAMSCALINGR1"):
                self.scaling_ranges.append( self.__parse_scaling(line, "param_scaling_range1", 2) )
                
            elif line.upper().startswith("PARAMSCALINGR2"):
                self.scaling_ranges.append( self.__parse_scaling(line, "param_scaling_range2", 2) )
                
            elif line.upper().startswith("PARAMSCALINGR3"):
                self.scaling_ranges.append( self.__parse_scaling(line, "param_scaling_range3", 2) )
                
            elif line.upper().startswith("KVOL"):
                self.kvol_given = True
                self.__parse_kvol(line)
                
            elif line.upper().startswith("THRESHOLD"):
                self.threshold_given = True
                self.__parse_positive_params(line, "threshold")
                
            elif line.upper().startswith("SIMBOX"):
                self.__parse_positive_params(line,'simbox', "float")
                
            elif line.upper().startswith("CRYSOL_PATH"):
                self.__parse_crysol_path(line)
                
            elif line.upper().startswith("GRIDRADIUS"):
                self.__parse_positive_params(line,'simboxradius', "float")
                
            elif line.upper().startswith("SAXSRADIUS"):
                self.__parse_positive_params(line,'saxsradius', "float")
            
            elif line.upper().startswith("RG_VAL"):
                self.__parse_positive_params(line,'rg_val', "float")
            
            elif line.upper().startswith("GRIDTYPE"):
                self.__parse_param_value(line, "simgridtype", ["cubic", "diamond"])
                
            elif line.startswith("SIMMETHOD"):
                self.__parse_param_value(line, "simmethod", ["SimulatedAnnealing", "Genetic", "ReplicaExchange"])
                
            elif line.upper().startswith("REPLICAEXCHANGE_FREQ"):
                self.__parse_positive_params(line, "replica_exchange_freq", "int")
                self.repl_exch_freq = True
                
            elif line.upper().startswith("MAXPOOLSIZE"):
                self.__parse_positive_params(line, "maxpoolsize", "int")
                
            elif line.upper().startswith("REPLICATEMPERATURES"):
                self.__parse_values_list(line, "replica_temps")
                
            elif line.upper().startswith("ANNTEMP"):
                self.__parse_positive_params(line,'anntemp', "float")
                
            elif line.upper().startswith("REDUCTMETHOD"):
                self.__parse_param_value(line, "reductmethod", ["Roulette", "Tournament", "Cutoff"])
                
            elif line.upper().startswith("MOVE_STATE"):
                self.__parse_movestate(line)
                
            elif line.upper().startswith("START_ORIENTATION"):
                self.__parse_bool_param(line, "start_orient")
                
            elif line.upper().startswith("IDENTIFY_DISORDERS"):
                self.__parse_bool_param(line, "identify_disorders")
                
            elif line.upper().startswith("COVALENT_BONDS"):
                self.__parse_covalent_links(line)
                

        fh.close()
        
        #parameter replicaexchangefreq not defined in config file
        if self.repl_exch_freq == False:
            self.replica_exchange_freq  = self.simul_steps/10
            if self.replica_exchange_freq == 0: self.replica_exchange_freq = 1

        self.__check_save_structnr(["struct_nr", "replica_exchange_freq"])
        
        self.__check_shape_descriptors()
        
        self.__check_mutual_clashes_options()
        
        self.__check_mut_freq_correctness()
        
        self.__set_outsave_mode()
                
        self.save_logfile()
                
    def __check_save_structnr(self, params):
        """
        """
        for param in params:
            if self.simul_steps == 0: pass
            elif getattr(self, param) > self.simul_steps:
                print getattr(self, param)
                raise ConfigError(str(param).upper()+" value cannot be larger than number of simulations!")
            
    def __check_mutual_clashes_options(self):
        if self.required_clashes_penalty == True and self.required_clashes_penalty_allatoms == True:
            raise ConfigError("Only one option can be applied CLASHES or CLASHES_ALLATOMS. \
                              Please change values into 0 for one of the parameters or comment one option in configuarion file")
        
    def __check_shape_descriptors(self):
        """
        """
        
        if self.shapedesc:
            if self.kvol_given == False and self.threshold_given == False and self.shapedesc == "map":
                raise ConfigError("You must provide kvol or threshold value when you provide density map as input")
        
        if self.kvol_given and self.threshold_given:
            raise ConfigError("Please provide only one of these two parameters: KVOL or THRESHOLD!")
        
        if self.threshold != None:
            self.kvol = None
        
        if self.is_freespace_defined == False:
            self.freespace_penalty = [0.,0.]
            
        if self.is_density_defined == False:
            self.density_penalty = [0.,0.]
        
        if self.density_penalty[0] != 0. and self.freespace_penalty[0] != 0.:
            print ConfigError("Scoring function will penalyze shape filling twice since you defined two parameters: DENSITY and MAP_FREESPACE!")
            logfile.write_file("Scoring function will penalyze shape filling twice since you defined two parameters: DENSITY and MAP_FREESPACE!\n")
            
        if self.simbox > 10:
            raise ConfigError("The size of the system you want to use is very large. Max simbox value is 10")
        
        #when no file with density map or ab initio model was provided but density filling or mapspace weights were provided
        if self.shapedesc == False and self.is_density_defined == True: #and (self.density_penalty[0] != 0):
            raise ConfigError("Map filling cannot be calculated when no shape descriptor was provided!\n")
            logfile.write_file("Map filling cannot be calculated when no shape descriptor was provided!")
        
        #print "@@@@@@", self.shapedesc, self.is_freespace_defined
        if self.shapedesc == False and self.is_freespace_defined == True: #(self.freespace_penalty[0] != 0):
            #print "****", self.shapedesc, self.is_freespace_defined
            raise ConfigError("Map filling cannot be calculated when no shape descriptor was provided!")
            logfile.write_file("Map filling cannot be calculated when no shape descriptor was provided!\n")
        
        if self.shapedesc == "map" or self.shapedesc == "saxs": self.shapedesc = True
        
        if self.is_chi2_defined == True or self.is_rg_penalty_defined == True:
            if not self.curve:
                raise ConfigError("To verify discrepancy with SAXS/SANS curves you must provide .dat file!")
            else:
                self.crysol_outfile = open("crysol_summary.txt", "w")
                
        if self.representation == "sphere" and self.restraints_penalty[0] == 0 and self.restraints_penalty[1] == 0:
            raise ConfigError ("To validate clashes between spheres PyRy3D calculates violation of distances betweeen\
                               spheres centres. To allow this option penalty for RESTRAINTS must be different than 0 0")
        
        if self.identify_disorders == False and self.simul_dd_freq != 0.:
            raise ConfigError ("You must allow PyRy3D to identify disorders to use simulate disorder mutation. Please set IDENTIFY_DISORDER parameter into True or disable simulate disorder mutation by setting SIMUL_DD_FREQ to 0 0 values")
    
        
    def __check_mut_freq_correctness(self):
        """
        """
        sumpar = float(self.rotation_freq + self.rotation_cov_freq + self.translation_freq + self.exchange_freq + self.exchangeandsample_freq\
                     + self.simul_dd_freq + self.translation_all_freq + self.rotation_all_freq \
                     + self.rotation_whole_freq)

        if sumpar == 0: raise ConfigError("Frequencies of mutations must sum up to 1. You provided %s"%(sumpar))
        if round(sumpar,1) > float(1.0):
            raise ConfigError("Frequencies of mutations must sum up to 1. You provided %s"%(sumpar))
        if round(sumpar,1) < float(1.0):
            self.rotation_freq        = self.rotation_freq/sumpar *1.0
            self.rotation_cov_freq    = self.rotation_cov_freq/sumpar *1.0
            self.translation_freq     = self.translation_freq/sumpar *1.0
            self.exchange_freq        = self.exchange_freq/sumpar *1.0
            self.exchangesample_freq  = self.exchangeandsample_freq/sumpar *1.0
            self.simul_dd_freq        = self.simul_dd_freq/sumpar *1.0
            self.translation_all_freq = self.translation_all_freq/sumpar *1.0
            self.rotation_all_freq    = self.rotation_all_freq/sumpar *1.0
            self.rotation_whole_freq  = self.rotation_whole_freq/sumpar *1.0
            
    def __set_outsave_mode(self):
        """
        """
        if self.simul_steps == 0: pass
        elif self.niter > self.simul_steps:
            raise ConfigError("Steps to write cannot be larger than number of simulation steps!")
        #self.out_steps = []
        
        if self.out_steps and self.niter and self.write_eachbetter:
            raise ConfigError("You can select only one of output save methods either WRITE_N_ITER or or WRITE_EACHBETTER or OUT_STEPS")
        
        if self.out_steps : self.save_res_mode = "outsteps"
        elif self.write_eachbetter   : self.save_res_mode = "eachbetter"
        elif self.niter   : self.save_res_mode = "niter"
        
        
    def save_logfile(self):
        """
        saves all alphabetically sorted config attributes to logfile
        """
        attribs = self.__dict__
        
        for a in sorted(attribs.keys()):
            logfile.write_file(str(a).upper()+"."*20+" "+str(attribs[a])+"\n")
            
    def set_movehistory_file(self, movehist_file):
        """
        """
        self.movehistory = movehist_file
        
class CovalentBond:
    def __init__(self, chains, at1=None, at2=None):
        self.atom1 = at1
        self.atom2 = at2
        self.chains = chains
        self.chains_indexes = []
        self.__check_atoms(self.atom1)
        self.__check_atoms(self.atom2)
        
    def __str__(self):
        return "%s %s %s" % (self.atom1, self.atom2, self.chains)
    
    def __check_atoms(self, atom):
        if len(atom) != 2: raise ConfigError("You haven't provided enough information about covalent bonds atoms")
        if not str(atom[0]).isdigit(): raise ConfigError ("Residue number should be number not string")
        if not str(atom[1]).isalpha(): raise ConfigError ("Atom name must be string, not number")
        
                
