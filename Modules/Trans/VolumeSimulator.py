#!/usr/bin/env python

from math                       import sin, cos, pi, sqrt, acos, atan2
from random                     import randint, uniform, sample, random
from numpy                      import array, dot
from numpy.random               import random
#import scipy.spatial

#PyRy3D
#from Modules.Trans.ComponentRepresentation import scipyconverter
from Modules.Input.PyRyStructure import PyryAtom
from Modules.Input.PyRyStructure import PyRyStructure
from Modules.Constans.Constans   import AMINOACIDS, NUCLEOTIDES
from Modules.Constans.Logfile    import logfile
from Modules.Error.Errors        import InputError
#biopython
import Bio
from Bio                        import PDB
from Bio.PDB                    import PDBParser, PDBIO, NeighborSearch
from Bio.PDB.Atom               import Atom
from Bio.PDB.Residue            import Residue
from Bio.PDB.Chain              import Chain
from Bio.PDB.Model              import Model
from Bio.PDB.Structure          import Structure


class point:
    x = 0.0
    y = 0.0
    z = 0.0
    
    def __init__(self, x=0, y=0, z=0):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        
    def __str__(self):
        return "[%f,%f,%f]" % (self.x, self.y, self.z)


def dist (p, q):
    return sqrt( ((p.x-q.x)**2) + ((p.y-q.y)**2) + ((p.z-q.z)**2) )
        
class Grapes (object):
    radius = 1.0
    masterRadius = 100.0
    #masterPoint
    #anchorPoint
    #grapes
    grapeCount = 10
    def __init__(self):
        self.masterPoint = point()
        self.generationMethod = "random" # available methods: "random", "chain", "anchored_chain"
        self.detectCollisions = None
        
    def __str__(self):
        return "%s %s %s %s %s" % (self.radius, self.masterPoint, self.masterRadius, self.grapeCount, self.generationMethod)


    def __set_anchor_coords(self, anchor_resi):
        """
        sets new coordinates of point object
        """
        anchor = point()
        anchor.x = anchor_resi[0]
        anchor.y = anchor_resi[1]
        anchor.z = anchor_resi[2]
        return anchor
        
    def __choose_random_anchor_position(self, anchor_residues):
        """
        sets position of first simulated atom position
        """
        anchors = []
        for anchor_resi in anchor_residues:
            anchor = self.__set_anchor_coords(anchor_resi)
            #randomly select translation axis
            axis = sample(["x", "y", "z"], 1)[0]
            #select translation direction
            vec = sample([self.radius*2, -self.radius*2], 1)[0]
            val = getattr(anchor, axis)
            setattr(anchor,axis, val + vec)
            anchors.append(anchor)
        return anchors

    def setRadius (self, r):
        self.radius = float(r);
        
    def setMasterRadius (self, r):
        self.masterRadius = float(r)

    def setMasterPoint(self, point, anchor=None):
        self.masterPoint = point
        if (not anchor is None):
            self.setAnchorPoint(anchor);
        
    def setAnchorPoint(self, point):
        self.anchorPoint = point
        
    def setGrapeCount (self, c):
        self.grapeCount = int(c)
        self.grapes = list([None]*c)
        
    def setMethod (self, method):
        self.generationMethod = method
        
    #   the argument should be a function that takes two arguments: center (point object), radius (double)
    def setCollisionDetectionFunction (self, func):
        self.detectCollisions = func
        
    def set_volume_simulation_parameters(self,fragment, component, struct, moltype, mass_centre):                                          
        """
        grapes_nr      - number of pseudoresidues to simulate
        sphere_radius  - radius of simulation sphere where pseudoresidues are created
        pa_radius      - radius of single pseudoresidue
        modeling_method- from selection of random and chain
        fragment_type  - from selection of terminal/internal/simulated_volume
        """
        self.setGrapeCount(len(fragment.sequence))
        self.setMasterRadius(fragment.max_sphere_radius)
        self.setRadius(fragment.radius)
        
        modeling_method = self.__choose_modeling_method(fragment.fragment_type, moltype)
        #print "Modeling method", modeling_method, fragment.fragment_type, pyrystruct.moltype
        self.setMethod(modeling_method)
        
        #collision detection only for components with some structure not for simulated_volumes
        if fragment.fragment_type != "simulated_volume":
            self.setCollisionDetectionFunction(self.check_if_collide)
              
            if struct: dd = list(struct.get_residues())
            else: dd = list(component.pyrystruct.struct.get_residues())
            dd.sort(key=lambda Residue: Residue.id[1])

#@TODO: when saving pdb files sort residues before!!

            if fragment.fragment_type == "nterm":
                anchor_resi = [self.__get_anchor_resi(component, struct, moltype, fragment, fragment.stop_pos+1)]

            elif fragment.fragment_type == "cterm":
                anchor_resi = [self.__get_anchor_resi(component, struct, moltype, fragment, fragment.start_pos-1)]

            elif fragment.fragment_type == "internal":
                anchor_resi  = [self.__get_anchor_resi(component, struct, moltype, fragment, fragment.start_pos-1),\
                                self.__get_anchor_resi(component, struct, moltype, fragment, fragment.stop_pos+1)]

            if len(anchor_resi) == 2:
                anchor_n = self.__choose_random_anchor_position(anchor_resi[0])
                anchor_c = self.__choose_random_anchor_position(anchor_resi[1])
                self.setMasterPoint(anchor_n[0], anchor_c[0])
            else:
                anchors = self.__choose_random_anchor_position(anchor_resi[0])
                self.setMasterPoint(anchors[0])
        
        else:
            masscenter = self.__set_anchor_coords(mass_centre)
            self.setMasterPoint(masscenter)
	    print "VOLUME SIMULATIONS IIIII", masscenter, mass_centre, self.masterPoint
                
    def __choose_modeling_method(self, fragment_type, moltype):
        """
        method decides on algorithm of fragment modeling dependable on
        molecule type (DNA, RNA, protein) and fragment location inside a structure
        (terminal, internal, simulated_volume)
        
        DNA in all cases is simulated as chain
        all internal fragments are chains
        terminal nucleic fragments as chains
        all protein fragments plus RNA simulated volume as random:
        """
        if fragment_type == "simulated_volume":
            if moltype == "dna": return "chain"
            else: return "random"
            
        elif fragment_type == "internal": return "anchored_chain"
        
        elif fragment_type == "cterm" or fragment_type == "nterm":
            if moltype == "protein": return "random"
            else: return "chain"        
        
#Funkcja losujaca punkt na SFERZE o srodku p i promieniu r
    def findSpherePoint (self, p, r) :
        alfa = uniform(0,1) * 2 * pi
        u = uniform (-1, 1)
        w = point()
        w.x = p.x + r * ( cos(alfa) * sqrt(1 - u**2) )
        w.y = p.y + r * ( sin(alfa) * sqrt(1 - u**2) )
        w.z = p.z + r * u
        return w

    def findBallPoint (self, p, r) :
        alfa = uniform(0,1) * pi
        beta = uniform(0,1) * 2.0 * pi
        q = uniform(0, r)
        w = point()
        w.x = p.x + q * sin(alfa) * cos(beta)
        w.y = p.y + q * sin(alfa) * sin(beta)
        w.z = p.z + q * cos(alfa)
        return w
    
    def findHemispherePoint (self, p, r, maxAngle, minAngle = 0.0, force=False) :
        fi = 2.0 * pi * uniform(0,1)
	if (force or maxAngle == minAngle):
	    teta = maxAngle
	else:
	    teta = acos (1 - uniform(0,1)*(1-cos(maxAngle)))
	    cnt = 0
	    while (maxAngle==pi and teta < minAngle and cnt < 10):
		teta = acos (1 - uniform(0,1)*(1-cos(maxAngle)))
		cnt += 1
	#fi and teta define a direction vector
        #print "angle chosen:", teta*180/pi
        w = point()
        w.x = p.x + r * cos(fi) * sin(teta)
        w.y = p.y + r * sin(fi) * sin(teta)
        w.z = p.z + r * cos(teta)
        return w
    
    def findChainPoint (self, p, r, maxAngle, minAngle = 0.0, force=False) :
	if (force or maxAngle == minAngle):
	    a = 0.0
	    teta = maxAngle
	else:
	    a = pi * 0.2 * uniform(-1, 1)
	    teta = maxAngle * uniform(-1, 1)
	    cnt = 0
	    while (teta < minAngle and cnt < 10):
		teta = maxAngle * uniform(-1, 1)
		cnt += 1
	w = point()
        w.x = p.x + r * sin(teta) * cos(a)
        w.y = p.y + r * sin(teta) * sin(a)
        w.z = p.z + r * cos(teta)
        return w
    
    def translate (self, v, history=True):
        for g in self.grapes :
            if (g is None):
                continue
            g.x += v[0]
            g.y += v[1]
            g.z += v[2]
        
        if (history):
            self.anchoredHistory.append (["T", v])
        
    def rotate (self, angle, axis, history=True):
        if   axis == 'X':  rot = array([[1,0,0],\
                                        [0, cos(angle), sin(angle)],\
                                        [0,-sin(angle), cos(angle)]])
        elif axis == 'Y':  rot = array([[cos(angle), 0, -sin(angle)],\
                                        [0,1,0], \
                                        [sin(angle), 0, cos(angle)]])  
        elif axis == 'Z':  rot = array([[cos(angle), sin(angle), 0],\
                                        [-sin(angle), cos(angle), 0],\
                                        [0,0,1]])
        for g in self.grapes:
            if (g is None):
                continue
            coords = [g.x, g.y, g.z]
            rotated_coordinates = dot(coords, rot) #multiply matrices
            g.x = rotated_coordinates[0]
            g.y = rotated_coordinates[1]
            g.z = rotated_coordinates[2]
            
        if (history):
            self.anchoredHistory.append(["R", angle, axis])

    def anchoredReposition(self, i) :
	if (i<1):
	    return
	#print "-----------------------"
	#print "repositioning ", i
        rootGrape = self.grapes[i-1]
	anchorGrape = self.grapes[i]
	#print "root", rootGrape, "anchor", anchorGrape

        v = [-rootGrape.x, -rootGrape.y, -rootGrape.z]
        self.translate(v)
	#print "translating by ", v
        
        xAngle = atan2(anchorGrape.y, anchorGrape.z)
	#print "rotating X by ", xAngle*180/pi
        self.rotate(xAngle, "X")
        yAngle = atan2(anchorGrape.x, anchorGrape.z)
	#print "rotating Y by ", yAngle*180/pi
        self.rotate(-yAngle, "Y")
	
	#print "root", rootGrape, "anchor", anchorGrape


    def checkChainDistance (self) :
        distance = dist(self.masterPoint, self.anchorPoint)
        return distance < ((self.grapeCount-1)*2.0*self.radius) # n-2 grapes need to fit into the chain plus two radiuses of the master and anchor
        
    def getMaxAngle (self, root, anchor, current) :
        distance = dist(root, anchor)
        anchorRadius = 2.0*self.radius*0.95
        rootRadius = (self.grapeCount-current)*2.0*self.radius
        if (distance > anchorRadius+rootRadius) :
            #print "max alph 0"
            return 0
        if (distance < abs(rootRadius-anchorRadius)) :
            #print "max alph 180"
            return pi
        a = (rootRadius**2 - anchorRadius**2 + distance**2 ) / (2.0*distance)
        height = sqrt(rootRadius**2 - a**2)
        maxAngle = atan2(height, a)
        if (distance < a):
            maxAngle = pi - maxAngle
        
        #print "current", current, "rootR", rootRadius, "dist", distance, "a", a, "h", height, "max alph:", maxAngle*180/pi
        return maxAngle
    
    def startAnchoredChain (self) :
        self.anchoredHistory = []
        self.grapes[1] = self.anchorPoint
        self.anchoredReposition(1)
	self.rotate(uniform(0, pi), "Z")
            
    def endAnchoredChain (self) :
        neg = lambda x: -x

        while (len(self.anchoredHistory)) :
            operation = self.anchoredHistory.pop()
            if (operation[0] == "T"):
                self.translate(map(neg, operation[1]), False)
            elif (operation[0] == "R"):
                self.rotate(-operation[1], operation[2], False)
    
    def collides (self, i):
        if (dist(self.masterPoint, self.grapes[i]) > (self.masterRadius - self.radius) ):
            return 1
        
        for j in range(0, i):
            if (dist(self.grapes[j], self.grapes[i]) < (self.radius + self.radius) ):
                return 1
        return 0
        
    def generate (self, component=None):
        """
        """
#        self.grapes[0] = self.findBallPoint(self.masterPoint, self.masterRadius-self.radius)
        self.grapes[0] = self.masterPoint
        if (self.generationMethod == "anchored_chain"):
            if (not self.checkChainDistance()):
                self.generationMethod = "chain"
                print "Warning! Anchored chain distance too big, switching to regular chain"
		logfile.write_file("Warning! Anchored chain distance too big, switching to regular chain\n")
            else:
                self.startAnchoredChain()
            
        for i in range(1, self.grapeCount):
            lIterCount = 0;
	    if (self.generationMethod == "anchored_chain"):
		self.anchoredReposition(i-1)

            while (lIterCount < 100):   #Try a hundred times before going ahead (last generated grape is left there, even if it collides)
                lIterCount += 1
                if (self.generationMethod == "random"):
                    n = randint(0, i-1) #attach to a random one
                elif self.generationMethod == "chain":
                    n = i-1 #attach new grape to the last one
		elif self.generationMethod == "anchored_chain":
                    #attention! i-2 can be equal -1
                    n = i-2 #attach new grape to the previous to last one
                    
                if (self.generationMethod == "anchored_chain"):
                    if (i == 1): #skip first step, we have the anchor point there
                        continue
		    angle = self.getMaxAngle(self.grapes[n], self.grapes[n+1], i)
                    self.grapes[i] = self.findChainPoint(self.grapes[n], 2*self.radius, 0.75*angle, 0.75*min(angle, pi/4), i==self.grapeCount-1)
                else:
                    self.grapes[i] = self.findSpherePoint(self.grapes[n], 2*self.radius)
                    
                if (self.collides(i) == 0):
                    if (self.detectCollisions != None):
                        if component:
                            if (self.detectCollisions(component, self.grapes[i], self.radius) == 0):
                                break
                    else:
                        break
	    #print "new grape", self.grapes[i]
        if (self.generationMethod == "anchored_chain"):
            self.endAnchoredChain()
        return self.grapes
    
    def clean_grapes(self):
        self.grapes = []
    
    def __get_anchor_resi(self, component,struct, moltype, fragment, *args): #res_number):
        """
        """
        anchor_residues = []
        if moltype == "protein":
            atom_name = "CA"
        else:
            atom_name = "C4'"
        for res_number in args:
            if struct:
                try: anchor_residues.append(list(struct.get_chains())[0][res_number][atom_name].coord)
                except: raise InputError("Residue %s in chain %s does not have CA atom"%(res_number, list(struct.get_chains())[0].id))
            else:
                try: anchor_residues.append(list(component.pyrystruct.struct.get_chains())[0][res_number][atom_name].coord)
                except: raise InputError("Residue %s in chain %s does not have CA atom"%(res_number, list(struct.get_chains())[0].id))
        return anchor_residues

    def check_if_collide(self, component, point, radius):
        """
        check whether aparticulat point does not collide with any components atoms
        
        Return all component's atoms that have at least one
        atom within radius of center for given point
        """
        ns               = NeighborSearch(list(component.get_atoms()))
        point_center     = array([point.x, point.y, point.z])
        #we assume that pseudoresidue radius ;lis 1.5A as a collistion detection area.
        found_collisions = len(ns.search(point_center, radius+1.5, "A"))
      
        return found_collisions

class Disordered_Fragment(object):
    
    def __init__(self, start_pos = None, stop_pos = None, sequence = None):
        if start_pos : self.start_pos = start_pos #residue number of first disordered residue
        else: self.start_pos   = 0
        if stop_pos  : self.stop_pos  = stop_pos  #residue number of last disordered residue
        else: self.stop_pos    = 0
        if sequence  : self.sequence  = sequence
        else: self.sequence    = ""              #sequence of disordered fragment
        self.radius            = 1.0              #pseudoatoms radius
        self.max_sphere_radius = 10.0             #radius of sphere defining volume simulation area
        self.pseudoresidues    = []               #list of residues objects
        self.fragment_type     = "internal"       #cterm/nterm/internal/simulated_volume
        self.fragment_lattice  = None             #lattice for disordered region structure (only residues, no atoms)
        
    def __str__(self):
        return "%s %s %s %s %s" % (self.start_pos, self.stop_pos, self.sequence, \
                                   self.radius, self.fragment_type)

    def add_component_structure(self, struct):
        """
        adds structure representing all component structure
        """
        self.structure = struct
        
    def add_fragment_structure(self, fragment):
        """
        adds piece of structure representing disordered region
        """
        self.fragment_lattice = fragment
        
    def add_pseudoatoms_to_structure(self, pseudoatoms, moltype):
        """
        """
        start_index = 0
        for pa in pseudoatoms:
        #    print "***", start_index -1, len(list(self.fragment_lattice.get_residues()))
            self.add_pa_to_structure(pa, list(self.fragment_lattice.get_residues())[start_index], moltype)
            start_index += 1
            
    def add_pa_to_structure(self, pa, resi, moltype):
        """
        """
        
        coord = array([pa.x, pa.y, pa.z]) #, "f")
        if moltype == "protein":
            new_atom = PyryAtom('CA', coord, 0, 1, ' ', ' CA', 1)
        else:
            new_atom = PyryAtom("C4'", coord, 0, 1, ' ', " C4'", 1)
            
        new_atom.assign_vdw()
        new_atom.assign_molweight()
        resi.add(new_atom)
        #print "Add PA to...", resi.id, new_atom.get_parent().id
        
#@TODO needs testing!!!!!
#@TODO need to renumber residues somehow!!

    def add_fragment_to_original_structure(self, component, structure, res_nr, fr_type):
        """
        normal - in direction from n to c term
        reverse - in direction from c to n term (for addition of residues on Nterm)
        """
        
#@TODO will have to be changed when multichain components come!!!
    #####################################    
        #chain = list(self.structure.get_chains())[0]
        if structure: chain = list(structure.get_chains())[0]
        else: chain = list(component.pyrystruct.struct.get_chains())[0]
        
        residues = list(self.fragment_lattice.get_residues())
        if fr_type == "nterm":
            residues.sort(key=lambda Residue: Residue.id[1]) #, reverse = True)
                
        for resi in residues:
            #print "WANNA ADD:  ", resi.id, res_nr, type
            resi.id = (" ", res_nr, " ")
            res_nr += 1
            chain.add(resi)
            self.add_pseudoresidue(resi)
        
    def add_pseudoresidue(self, pr):
        """
        adds pseudoresidue object to list of pseudoresidues
        """
        self.pseudoresidues.append(pr)
        
    def build_structure(self, sequence, start_index, moltype):
        """
           builds new structure composed of atom_name atoms only
           for nucleic acids these are C4', for proteins CA
        """
        new_chain = list(self.fragment_lattice.get_chains())[0]
        for resi in sequence:
            #print "resi", resi, start_index
            resi_id = (" ", start_index, " ")
            #resi_name = resi
###############33
	    if moltype == "protein":
		resi_name = AMINOACIDS[resi.upper()]
	    else:
		resi_name = NUCLEOTIDES[resi.upper()]
###############33
            new_resi = Residue(resi_id, resi_name, " ")
            new_chain.add(new_resi)
            start_index += 1
            
    def clean_fragment(self):
        self.fragment_lattice = None
        
    def remove_pseudoresidues(self, structure):
        """
        removes pseudoresidues simulated during previous mutation in order
        to prepare conditions for new simulated pseudoresidues to be attached and
        scored
        """
#############################
##self.structure into structure
#@TODO must be changed when hybrids components will be considered
        remove = []
        if self.fragment_type != "simulated_volume":
            if structure:
                chain = list(structure.get_chains())[0]
                for resi in chain:
                    for r in self.pseudoresidues:
                        if r.id[1] == resi.id[1]:
                            remove.append(resi)
                            #print "wanna detach", resi.id, len(list(self.structure.get_residues()))
                            break
            for resi in remove:
                chain.detach_child(resi.id)
        
        self.pseudoresidues = []        
        #self.fragment_lattice = None

    def create_new_chain(self, id):
        """
        """
        self.fragment_lattice = Structure(id) 
        my_model = Model(0)
        self.fragment_lattice.add(my_model)
        my_chain = Chain(id)
        my_model.add(my_chain) #what if more chains in one component?
        
    def create_simulated_volume(self, start_pos, fasta_seq, struct):
        """
        method to create new instance of Disordered_Fragment class to represent
        regions with not assigned atom coordinates
        """
        self.set_fragment_type("simulated_volume")
        self.create_new_chain(struct.chain)
        self.build_structure(fasta_seq, start_pos, struct.moltype)
        self.set_fragment_sequence(fasta_seq)
        self.get_pseudoatom_radius(struct)
##@TODO-CHECK: how to assess radius of simulation sphere??
        self.calculate_max_sphere_radius(struct)


############################################33
#@TODO-CHECK: wouldn't one fragment type be enough??

    def create_simulated_fragment(self, struct, fasta_seq):
        """
        method to create set attributes of disordered fragment instance
        """
#@TODO-CHECK: calculate max sphere radius according to moltype and resi number!!
        self.get_pseudoatom_radius(struct)      
        self.__check_fragment_type(fasta_seq)
        self.calculate_max_sphere_radius(struct)

        
    def __check_fragment_type(self, fasta_seq):
        """
        from selection of:
           terminal
           internal
        """
	################333
	#if len(fasta_seq) >30:
	#fragments longer than 30 residues are simulated as grapes.
	#################33
        if self.stop_pos == len(fasta_seq): self.fragment_type = "cterm"
        elif self.stop_pos < len(fasta_seq)-1 and self.start_pos > 1:
            self.fragment_type = "internal"
	    if self.stop_pos - self.start_pos >= 30:
	        self.fragment_type = "cterm"
	        print "internal fragment simulated as GRAPE", self.start_pos, self.stop_pos
        elif self.start_pos == 1: self.fragment_type = "nterm"
	
            
###########################################
        
    def get_pseudoatom_radius(self, struct):
        """
        radius is averaged distance between CA or P atoms in ideal helix/2 to
        represent volume of single CA/P atom
        """
        if   struct.moltype.lower() == "protein": self.radius = 1.9 #or 1.72 as C atom
        elif struct.moltype.upper() == "DNA"    : self.radius = 3.6
        elif struct.moltype.upper() == "RNA"    : self.radius = 3.0
                    
    def calculate_max_sphere_radius(self, struct):
        """
           returns average residue radius for a particular component type (in Angstrooms)
           
           radius given (3.8, mean 6.8, 5.8) is averaged distance between CA or C4'
           atoms in ideal helix
        """
        if   struct.moltype.lower() == "protein":
            self.__set_max_sphere_radius_for_moltype(3.8)
        
        elif struct.moltype.upper() == "DNA"    :
            self.__set_max_sphere_radius_for_moltype(7.0)
         
        elif struct.moltype.upper() == "RNA"    :
            self.__set_max_sphere_radius_for_moltype(6.0)
                
    def __set_max_sphere_radius_for_moltype(self, mol_radius):
        """
        """
        if self.fragment_type   == "cterm" or  self.fragment_type   == "nterm":
            self.max_sphere_radius = (len(self.sequence)*mol_radius)/2
        elif self.fragment_type    == "simulated_volume":
            self.max_sphere_radius = (len(self.sequence)*mol_radius)
            
    def set_max_sphere_radius(self, area_radius):
        """
        sets radius of simulation area for Volume Simulator
        """
        self.max_sphere_radius = area_radius

    def set_anchor_residues(self, resi1, resi2 = None):
        """
        defines first residue for volume simulator - here simulation should start
        if given defines last residue for volume simulator - here simulation should finish
        """
        self.start_resi = resi1
        self.end_resi   = resi2
        
    def set_modeling_disordered_fragment(self, component, struct, chain):
        """
#@TODO: podzial na 2 osobne funkcje: modelowanie fragmentow i modelowanie objetosci
#@TODO: uporzadkowac simulate fragments!! dodac symulacje dla fragmentow srodkowych i nkoncowych
#@TODO: sprawdzic numeracje dodawanych przez symulator reszt
        """
        if self.fragment_type == "simulated_volume":
#@TODO: must be changed into nicer way!!
            if struct:
                self.remove_pseudoresidues(struct) # pyrystruct.struct#remove old Pseudoatoms positions
            else:
                self.remove_pseudoresidues(component.pyrystruct.struct)
            self.create_new_chain(chain) #pyrystruct.chain
            self.build_structure(self.sequence, self.start_pos, component.pyrystruct.moltype) 
        else:
            if struct:
                self.remove_pseudoresidues(struct)
            else:
                self.remove_pseudoresidues(component.pyrystruct.struct)
            self.create_new_chain(chain)
            self.build_structure(self.sequence, self.start_pos, component.pyrystruct.moltype)
            #self.add_component_structure(pyrystruct.struct)
            
    def set_fragment_sequence(self, seq):
        """
        sequence of disordered fragment
        """
        self.sequence = seq

    def set_pseudoatom_radius(self, radius):
        """
        sets pseudoatom radius; different for nucleotide and amino acids; in Angstroms
        """
        self.radius = radius

    def set_fragment_type(self, frag_type):
        """
        sets fragment type. Can be internal when disordered region is inside the component's structure
        or terminal when it is located on C or N termini
        """
        self.fragment_type = frag_type


#@TODO: helping function: for tests isues only; to be removed!!
def save_pdb(struct, name):
    print "test",len(list(struct.get_residues()))
    for resi in struct.get_residues():
        print resi.id, resi.resname
    out = PDBIO()
    out.set_structure(struct)
    out.save(str(name)+'volume_simulator.pdb')


if __name__=='__main__':

    p = PDBParser(PERMISSIVE=False, QUIET=True)
    st = p.get_structure("Zfull.pdb", "Zfull.pdb")
    component = Component()
    component.pyrystruct = st
    
    
    dd_frag = DisorderedFragment()
    dd_frag.set_modeling_disordered_fragment(component.pyrystruct)
    mass_centre = [0,0,0]
    
    g = Grapes()
    g.set_volume_simulation_parameters(dd_frag, component.pyrystruct, mass_centre)
    
    res = g.generate()
    
    dd_frag.add_pseudoatoms_to_structure(res, component.pyrystruct.moltype)
    
    dd_frag.add_fragment_to_original_structure(94, "chain")   #(start_pos)

    
    last_atom = list(st.get_atoms())[-1].coord
    print "last atom", last_atom
    
    la = point()
    la.x = last_atom[0]
    la.y = last_atom[1]
    la.z = last_atom[2]
    
    
    g = Grapes()
    g.setGrapeCount(10)
    g.setMasterRadius(30.0)
    g.setRadius(1.0)
    g.set_fragment(frag)
    g.setMethod("anchored_chain")
    g.component = st
    g.setCollisionDetectionFunction(g.check_if_collide)
#       You can add a collision detection function here if You have a Component object, with it's find_collisions function
#       see Component.py line 416 for an example
    g.setMasterPoint(la)
    res = g.generate()
    for grape in res:
        print "SPHERE, ", grape.x, ", ", grape.y, ", ", grape.z,", ", 1.0, ","

    #frag.add_pseudoatoms_to_structure(res, 0)
    frag.add_pseudoatoms_to_pseudoresidues(res)
    save_pdb(frag.structure, "bla1111")
    #frag.remove_old_pseudoatoms(5,14)
    #save_pdb(frag.structure, "bla2222")
    
    #comp  = Component()
    #ps = PyRyStructure(frag.structure)
    #comp.pyrystruct = ps
    #comp.set_center()
    #comp.rotate(30, "X")
    #save_pdb(comp.pyrystruct.struct, "dddd")

#przed symulacja usun wszystkie pseudoreszty (detach_child) - poza atomami fixed!!  --done
#to klasa fragmnet musi decydowac jakie regiony maja byc symuowne i w jakim trybie
#kazdy punkt jako pseudoreszta musi byc dodany do struktury komponentu!/dodaj reszte zbudowana z jednego atomu o okreslonym centrum



