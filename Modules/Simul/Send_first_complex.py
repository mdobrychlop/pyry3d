
from Modules.Simul.cpp				import pyry3d_cpp
from sys				import exit

def send_coords(first_complex, complex_index):
	
	"""
		send first_complex to cpp module,
		pdb - pdb format!!
    """	
	pdb = '\n'
	at_nr = 0
	for component in first_complex.components:
		for at in component.pyrystruct.struct.get_atoms():
			at_nr += 1
			residue = at.get_parent()
			chain = residue.get_parent()			
			pdb += "ATOM".ljust(6,' ') + repr(at_nr).rjust(5,' ') + " " +\
			str(at.fullname).ljust(4,' ') + " " + (str(residue.resname)) +\
			" " + (str(chain.id))
			pdb += (str(residue.id[1])).rjust(4, ' ') + "    " +\
			'{0:.3f}'.format(round(at.coord[0],3)).rjust(8) +\
			'{0:.3f}'.format(round(at.coord[1],3)).rjust(8) +\
			'{0:.3f}'.format(round(at.coord[2],3)).rjust(8)+\
			'{0:.1f}'.format(at.occupancy).rjust(5) +\
			" " +'{0:.2f}'.format(at.bfactor).rjust(5) + "          "\
			+ str(at.element).rjust(2) + " \n"
	pyry3d_cpp.parse_pdb(pdb)
	
def send_first_grid(mDensityMap):
	"""
		send whole Grid data to cpp module,
		more about Grid and grid_cells in complex.hpp
    """	
	grid = mDensityMap.mapgrid
	simulbox = mDensityMap.simulbox
	#~ print "grid.mapcells: " + str(len(grid.mapcells))
	#~ print "sim.mDensityMap..density_sum: " + str(sim.mDensityMap.density_sum)
	grid_cells_index = []
	grid_cells_ismap = []
	grid_cells_density = []
	grid_cells_x = []
	grid_cells_y = []
	grid_cells_z = []
	for i in grid.grid_cells:
		grid_cells_index.append(i.index)
		grid_cells_ismap.append(i.ismap)
		grid_cells_density.append(i.density)
		grid_cells_x.append(i.coord[0])
		grid_cells_y.append(i.coord[1])
		grid_cells_z.append(i.coord[2])

	pyry3d_cpp.send_grid(len(grid.mapcells), mDensityMap.density_sum,
	simulbox.xmin, simulbox.xmax, simulbox.ymin, simulbox.ymax,\
	simulbox.zmin, simulbox.zmax, grid.xmin, grid.xmax, grid.ymin, grid.ymax,\
	grid.zmin, grid.zmax, grid.radius, grid.overlap, grid_cells_index,\
	grid_cells_density, grid_cells_ismap, grid_cells_x, grid_cells_y, grid_cells_z)

def send_first_complex(first_complex, complex_index):
	"""
		send first_complex data to cpp module,
		ready for mutations
		(with statistics calculated in python)	
		more about Complex and Components in complex.hpp
    """	
	taken_mapcells = len(first_complex.taken_mapcells)
	pyry3d_cpp.set_changed_component(complex_index, first_complex.changes)
	# uwaga rzutowanie z float na int (jawne)
	pyry3d_cpp.set_complex_data(complex_index, first_complex.alfa_atoms_nr,\
	first_complex.all_clashes, int(first_complex.all_filled_cells),\
	int(first_complex.all_outbox_atoms), 0,\
	first_complex.clashes, first_complex.density, first_complex.is_complex_outbox,\
	first_complex.is_complex_outmap, first_complex.mapcells_nr,\
	first_complex.pairs_atoms, first_complex.restraints, first_complex.simulation_score,\
	first_complex.symmetry, first_complex.taken_densities, taken_mapcells)
	
	is_outbox = []
	is_outmap = []
	mapcells_nr = []
	outbox_atoms = []
	for component in first_complex.components:
		is_outbox.append(component.is_outbox)
		is_outmap.append(component.is_outmap)
		mapcells_nr.append( int(component.mapcells_nr) )
		outbox_atoms.append(component.outbox_atoms)

	pyry3d_cpp.set_components_data(complex_index, is_outbox,\
	is_outmap, mapcells_nr, outbox_atoms)

def send_first_clashes(first_complex, complex_index):
	"""
		send to cpp module information about clashes calculated in python
    """
	clashes_nr = []
	for i in first_complex.pairs:
		clashes_nr.append(i.pair[0])
		clashes_nr.append(i.pair[1])
		clashes_nr.append(i.clashes_nr)
	
	pyry3d_cpp.set_clashes_nr(complex_index, clashes_nr)
