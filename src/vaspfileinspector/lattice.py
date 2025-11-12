# -*- coding: utf-8 -*-
import numpy as np
from math import (sin,cos,sqrt)
import spglib
import sys
from vaspfileinspector.common import (R2D,D2R)

class Lattice:

	def __init__(self,H):

		self.a = 0
		self.b = 0
		self.c = 0
		self.alpha = 0
		self.beta = 0
		self.gamma = 0

		self.volume = 0

		self.symIntlSymb = "P1"
		self.spgNumber = 1
		self.bravais = "unknown"

		self.H = H

		self.alat = 1.0
		self.a1 = [float(x) for x in H[0]]
		self.a2 = [float(x) for x in H[1]]
		self.a3 = [float(x) for x in H[2]]

		self.box2cell()
		self.set_volume()

		## finish construction

	def analyze_symmetry(self,cell,sp):

		dataset = spglib.get_symmetry_dataset(cell,sp)

		self.spgNumber = dataset['number']

		self.symIntlSymb = dataset['international']

		if self.spgNumber <= 2:
			self.bravais = "triclinic"
		elif self.spgNumber <= 15:
			self.bravais = "monoclinic"
		elif self.spgNumber <= 74:
			self.bravais = "orthorhombic"
		elif self.spgNumber <= 142:
			self.bravais = "tetragonal"
		elif self.spgNumber <= 167:
			self.bravais = "trigonal"
		elif self.spgNumber <= 194:
			self.bravais = "hexagonal"
		elif self.spgNumber <= 230:
			self.bravais = "cubic"
		else:
			self.bravais = "unknown"            

	def show_info(self,parameters,cell,sp):

		if parameters.save:
			name = parameters.compound + ".cell"
			out = open(name,'w')
		else:
			out = sys.stdout
		try:
			
			if out is sys.stdout:
				out.write('\n' + "/*-- Cell --*/" + '\n')
			else:
				out.write("/*-- Cell --*/" + '\n')
			out.write("Structure     = %s    " % parameters.FILE + '\n')
			out.write("Compound      = %s    " % parameters.compound + '\n')
			out.write("symprec       = %3.0e " % sp + '\n\n')
			# out.write("** Cell (H) matrix:" + '\n')

			out.write("        | %-5f %5f %5f |" % (self.a1[0],self.a1[1],self.a1[2]) + '\n')
			out.write("  H  =  | %-5f %5f %5f |" % (self.a2[0],self.a2[1],self.a2[2]) + '\n')
			out.write("        | %-5f %5f %5f |" % (self.a3[0],self.a3[1],self.a3[2]) + '\n\n')
			# out.write("** Cell vectors and angles:" + '\n')
			out.write("Volume   = %-10f (Å^3)" % self.get_volume() + '\n')
			out.write("a        = %-10f (Å)" % self.a + '\n')
			out.write("b        = %-10f (Å)" % self.b + '\n')
			out.write("c        = %-10f (Å)" % self.c + '\n')
			out.write("alpha    = %-10f (°)" % (self.alpha*R2D) + '\n')
			out.write("beta     = %-10f (°)" % (self.beta*R2D) + '\n')
			out.write("gamma    = %-10f (°)" % (self.gamma*R2D) + '\n')
			out.write("symmetry = %-s %-i  " % (self.symIntlSymb,self.spgNumber) + '\n')
			out.write("bravais  = %-s" % self.bravais + '\n')
		finally:
			if parameters.save:
				out.close()


	def scale(self,a):
		self.alat = a
		self.a1 = [a*x for x in self.a1]
		self.a2 = [a*x for x in self.a2]
		self.a3 = [a*x for x in self.a3]


	def get_volume(self):
		if self.volume is None:
			self.set_volume()
		else:
			return self.volume

	def set_volume(self):
		ca = cos(self.alpha);sa = sin(self.alpha)
		cb = cos(self.beta); sb = sin(self.beta)
		cg = cos(self.gamma);sg = sin(self.gamma)

		a = self.a
		b = self.b
		c = self.c

		self.volume = a*b*c * sqrt(1.0- (ca*ca) - (cb*cb) - (cg*cg) + 2*(ca * cb * cg))

	# Convert lower triangular matrix "H",
	#      |  H[0][0]     0        0      |
	# H =  |  H[1][0]  H[1][1]     0      |
	#      |  H[2][0]  H[2][1]  H[2][2]   |
	#     into conventional lattice vectors defined by 
	#     "a" "b" "c" "alpha" "beta" "gamma"
	def box2cell(self):
		self.a = np.linalg.norm( self.H[0] )
		self.b = np.linalg.norm( self.H[1] )
		self.c = np.linalg.norm( self.H[2] )

		self.alpha = self._get_angle( self.H[1],self.H[2])
		self.beta = self._get_angle( self.H[0],self.H[2])
		self.gamma = self._get_angle( self.H[0],self.H[1])


	def _get_angle(self,u,v):
		c = np.dot(u,v)/np.linalg.norm(u)/np.linalg.norm(v) # -> cosine of the angle
		angle = np.arccos(np.clip(c, -1, 1)) # if you really want the angle
		return angle
