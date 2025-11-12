# -*- coding: utf-8 -*-

import numpy as np
import math
from vaspfileinspector.common import Point
import sys


class Neighbors:
	def __init__(self,rcut=0):

		self.indicies = []
		self.neighbors = []
		self.bonds = []
		
		if rcut == 0:
			self.search = True
			self.rcut = 0
		else:
			self.search = False
			self.rcut = rcut

		self.nbonds =0

		self.minPair = None

		self.minBond = 10


	def show_info( self,parameters,atoms,depth=1 ):

		bondList = self.get_bond_list()

		indexer,neighborList = self.get_nn_list()

		minPair,minBond = self.get_min_pair()

		if parameters.save:
			name = parameters.compound + ".bonds"
			out = open(name,'w')
		else:
			out = sys.stdout
		try:
			if depth > 0:
				if out is sys.stdout:
					out.write('\n' + "/*-- Bonds --*/" + '\n')
				else:
					out.write("/*-- Bonds --*/" + '\n')
				out.write("Structure     = %s    " % parameters.FILE + '\n')
				out.write("Compound      = %s    " % atoms.get_compound() + '\n')
				out.write("NBonds        = %i    " % self.nbonds + '\n')
				out.write("Search Radius = %f  (Å)  " % self.rcut + '\n')
				out.write("Minimum Bond  = %f  (Å)  " % minBond + '\n')
				out.write("Minimum Pair  = %s    " % minPair + '\n')

			if depth > 1:
				out.write("/*-- Neighbors --*/" + '\n')
				for i in range(len(atoms.x)):
					nneighs = indexer[i+1] - indexer[i]
					out.write(" %s%i atom(#%i) has %i neighbors:" % (atoms.symbols[i],atoms.ids[i],i+1,nneighs) + '\n')
					for pj in range(indexer[i],indexer[i+1]):
						if i != pj:
							out.write("   %s%i-%s%i = %5f" % (atoms.symbols[i],atoms.ids[i],atoms.symbols[neighborList[pj]],atoms.ids[neighborList[pj]],bondList[pj]) + '\n')
		finally:
			if parameters.save:
				out.close()

	def set_min_pair(self,bond,i,s1,j,s2):
		last = self.minBond
		self.minBond = min(self.minBond,bond)

		if self.minBond < last:
			self.minPair = s1 + str(i) + "-" +  s2 + str(j)

	def get_min_pair(self):
		return self.minPair,self.minBond

	def get_bond_list(self):
		return self.bonds

	def get_nn_list(self):
		return self.indicies,self.neighbors

	def distance( self,atomi,atomj ):
		d = math.sqrt(math.pow(atomi[0]-atomj[0],2) + math.pow(atomi[1]-atomj[1],2) + math.pow(atomi[2]-atomj[2],2))
		return d

	def find(self,atoms,lattice,species,rcut):
		found = self.build_list( atoms,lattice,species,rcut )
		if self.search:
			while not found:
				found = self.build_list( atoms,lattice,species,rcut )
				if not found:
					rcut += 0.2
				else:
					break

	def species_index( self,species,j ):

		count = 1
		symbol = species[j]
		for i in range(len(species)):
			if symbol == species[i]:
				if i == j:
					break
				else:
					count += 1
		return count


	def build_list(self,atoms,lattice,species,rcut):
		# atoms[id][x,y,z]
		# lattice[ai][x,y,z]

		idx = [-1,0,1]
		found = False

		self.rcut = rcut

		nsize = 200*len(atoms)

		for i in range(len(atoms)+1):
			self.indicies.append(0)
		for i in range(10*len(atoms)):
			self.neighbors.append(0)
			self.bonds.append(0)

		atomj = []
		atomj.append(0)
		atomj.append(0)
		atomj.append(0)

		for iz in idx:
			for iy in idx:
				for ix in idx:
					shiftX = (ix * lattice.a1[0]) + (iy * lattice.a2[0]) + (iz * lattice.a3[0])
					shiftY = (ix * lattice.a1[1]) + (iy * lattice.a2[1]) + (iz * lattice.a3[1])
					shiftZ = (ix * lattice.a1[2]) + (iy * lattice.a2[2]) + (iz * lattice.a3[2])

					pointer = 0
					for i in range(len(atoms)):
						for j in range(len(atoms)):
							atomj[0] = atoms[j][0] + shiftX
							atomj[1] = atoms[j][1] + shiftY
							atomj[2] = atoms[j][2] + shiftZ

							bij = self.distance( atoms[i],atoms[j] )

							if bij != 0 and bij <= rcut:
								found = True
								jc =0
								if ix == 0 and iy == 0 and iz == 0:
									self.nbonds += 1
								jc = self.species_index( species,j )
								ic = self.species_index( species,i )

								self.set_min_pair(bij,ic,species[i],jc,species[j])
								self.neighbors[pointer] = j
								self.bonds[pointer] = bij
								pointer += 1
								if pointer >= nsize:
									print("warning: neighbor overflow")
									nsize += nsize/5
									for i in range(nsize/5):
										self.neighbors.append(0)
						self.indicies[i+1] = pointer

		return found

