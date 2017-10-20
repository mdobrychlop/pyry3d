#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# 
# www.genesilico.pl 
#

from Modules.Simul.cpp				        import pyry3d_cpp

def cpp_cutoff(sim):
	while( pyry3d_cpp.simul_size() >  sim.mMaxPoolSize):
		pyry3d_cpp.pop()
	pyry3d_cpp.set_pool_size( sim.mMaxPoolSize )



def cpp_roulette(sim):
	pyry3d_cpp.roulette(sim.mMaxTemperature, sim.mMaxPoolSize)


def cpp_tournament(sim):
	pyry3d_cpp.tournament(sim.mMaxPoolSize)
	