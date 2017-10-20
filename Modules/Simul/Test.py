from Bio.PDB.Structure              import Structure
from Bio.PDB.Model                  import Model
from Bio.PDB                        import PDBIO
from Modules.Simul.cpp				import pyry3d_cpp


def complex_save(given_complex, i, path):
	
	s = Structure(i) 
	my_model = Model(0)
	s.add(my_model)
	for component in given_complex.components:
		my_model.add(component.pyrystruct.struct[0][component.pyrystruct.chain])
	out = PDBIO()
	out.set_structure(s)
	out.save(path)
	return path	

def write_test_data(first_complex, index):
	
	mapcells_nr = []
	for i in first_complex.components:
		mapcells_nr.append(i.mapcells_nr)

	complex_data = []
	complex_data.append(first_complex.alfa_atoms_nr)
	complex_data.append(first_complex.all_clashes)
	complex_data.append(first_complex.all_filled_cells)
	complex_data.append(int(first_complex.all_outbox_atoms))
	complex_data.append(first_complex.all_restraints)
	complex_data.append(first_complex.clashes)
	complex_data.append(first_complex.density)

	if first_complex.is_complex_outbox:
		complex_data.append(1)
	else:
		complex_data.append(0)
	if first_complex.is_complex_outmap:
		complex_data.append(1)
	else:
		complex_data.append(0)

	complex_data.append(first_complex.mapcells_nr)	
	complex_data.append(first_complex.pairs_atoms)
	complex_data.append(first_complex.restraints)
	complex_data.append(first_complex.simulation_score)
	complex_data.append(first_complex.symmetry)
	complex_data.append(first_complex.taken_densities)
	complex_data.append(len(first_complex.taken_mapcells))

	components_data = []
	taken_mapcells = []
	
	for components in first_complex.components:
		taken_mapcells.extend(components.taken_mapcells)
		taken_mapcells.append(-1)
		components_data.append(len(components.alfa_atoms))

		if components.is_outbox:
			components_data.append(1)
		else:
			components_data.append(0)

		if components.is_outmap:
			components_data.append(1)
		else:
			components_data.append(0)
			
		components_data.append(components.mapcells_nr)	
		components_data.append(components.outbox_atoms)
		components_data.append(components.taken_density)


	pyry3d_cpp.write_get_complex_data_float(complex_data)
	pyry3d_cpp.write_get_components_data_float(components_data)
	complex_save(first_complex, index, "test_coords")
	pyry3d_cpp.write_changed("test_changed", first_complex.changes)
	pyry3d_cpp.write_mapcells_nr("test_mapcells_nr", mapcells_nr)
	pyry3d_cpp.write_taken_mapcells("test_taken_mapcells", taken_mapcells)

	clashes_nr = []			
	for i in first_complex.pairs:
		clashes_nr.append(i.pair[0])
		clashes_nr.append(i.pair[1])
		clashes_nr.append(i.clashes_nr)
	#~ print "tuz przed mutacja clashes: " + str(i.pair[0]) + " " + str(i.pair[1]) + " " + str(i.clashes_nr)
	pyry3d_cpp.write_clashes_nr_components(clashes_nr)
	
	
def read_test_data(index):
	pyry3d_cpp.test_get_complex_data(index)
	pyry3d_cpp.test_get_components_data(index)
	pyry3d_cpp.test_get_coords("test_coords", index)
	pyry3d_cpp.test_get_changed_component("test_changed", index)
	pyry3d_cpp.test_set_taken_mapcells(index)#from file
	pyry3d_cpp.test_set_clashes_nr(index)
	
def send_taken_mapcells(lNewComplex):
	res = []
	for i in lNewComplex.changes:
		component = lNewComplex.components[i]
		res.append(len(component.taken_mapcells))
	pyry3d_cpp.send_taken_mapcells(res)
