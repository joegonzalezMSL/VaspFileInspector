# -*- coding: utf-8 -*-

import sys
import argparse
import numpy as np
import textwrap

# import reader
# import common
# from neighbors import *
# from lattice import *
# from atoms import *
from vaspfileinspector import common, reader
from vaspfileinspector.neighbors import *
from vaspfileinspector.lattice import *
from vaspfileinspector.atoms import *

import spglib

## Author: Joe Gonzalez, Department of Physics, University of South Florida
## Version 3: 08/2017

def get_arguments(argv):

	cli = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''\
           Utility to analyze POSCAR/CONTCAR structure file
           By default print all information; unit cell, bonds, atoms ...
           Selectively print info by using optional arguments.
         -----------------------------------------------------------------
         '''),
    epilog=textwrap.dedent('''\
         examples:
         -----------------------------------------------------------------
            %(prog)s POSCAR > structure-data.nfo
            %(prog)s -p POSCAR
            %(prog)s -f POSCAR -bc > bonding-data.nfo 
            %(prog)s --debug --radius=2 --neighbors --save sns2.vasp --primitive
            %(prog)s -vvvvvnacrs 3 --tolerance=1e-3 mos2.contcar 
            %(prog)s -vvnacrs 3 mos2.contcar -t 0.01
            %(prog)s -gnacrs3 mos2.contcar -t0.1
         '''))

	cli.add_argument("FILE",help="input file containig data to process",type=str)
	cli.add_argument("-a","--atoms",dest="printAtoms",help="print atomic scale info. ntotal, types ...",action="store_true")
	cli.add_argument("-b","--bonds",dest="printBonds",help="print bonding info. total, species connectivity",action="store_true")
	cli.add_argument("-c","--cell",dest="printCell",help="print info on the unit cell; a,b,c volume ... ",action="store_true")
	cli.add_argument("-n","--neighbors",dest="printNlist",help="print bonding info and nearest neighbor information. n^2 recursive search for nearest neighbors stop when at least (1) bond is made, use this for assigning all neighbors.",action="store_true")
	cli.add_argument("-p","--primitive", dest="getPrimitive",help="if possible,reduce convetional cell to primitive unit cell",action="store_true")
	cli.add_argument("-r","--radius=", dest="rcut",help="search radius for considering an atom as nearest neighbor,(default = %(default)s Å)",default=0.0,type=float)
	cli.add_argument("-s","--save",dest="save",help="save the computed data to a file <stoich>.[bonds,atoms,cell]. ex: P1S2H2 -> P1S2H2.bonds P1S2H2.atoms, default == do not save, print to stdout",action="store_true")
	cli.add_argument("-t","--tolerance=", dest="symprec",help="precision in determining crystal symmetry in Cartesian coordinates,(default = %(default)s Å)",default=0.05,type=float)
	cli.add_argument("-v", dest="verb",help="increase output verbosity",default=0,action="count")
	cli.add_argument("--debug", dest="verb",help="extensive info, equivalent to \"-vvvv\"",action="store_true")
	cli.add_argument('--version', action='version', version='%(prog)s 4.0.')
	args = cli.parse_args()
	# if args.verb>0: print "verbosity level: ", args.verb
	# if args.verb>1: print "verbosity level: ", args.verb

	return args


def main():

	argv = sys.argv[1:]
	parameters = get_arguments( argv )

	if len(argv) == 1:
		parameters.printNlist = True
		parameters.printBonds = True
		parameters.printAtoms = True
		parameters.printCell = True

	# data[0] -> lattice
	# data[1] -> positions  in angstrom
	# data[2] -> species
	# data[3] -> ntypes
	# data[4] -> (a)tomic (n)umber(s)
	data = reader.read_vasp( parameters.FILE )

	# create instance of Neighbors
	nn = Neighbors(parameters.rcut)

	# store data read from structure file
	lattice = Lattice(data[0])
	positions = data[1]
	species = data[2]
	ntypes = data[3]
	ans = data[4]

	# read in the data and store in the "atom" object
	atoms = Atoms(lattice.H,lattice.volume,species,positions,ans,ntypes,fractional=False)
	parameters.compound = atoms.get_compound()

	# combine structural data, used for spacegroup/symmetry info
	cell = (lattice.H,atoms.xs,ans)

	# determine the symmetry, spacegroup, bravais, etc. 
	# of the lattice
	lattice.analyze_symmetry(cell,parameters.symprec)

	# if we want bond level info, build the neighbor list first
	if parameters.printBonds or parameters.printNlist:
		nn.find(atoms.x,lattice,species,parameters.rcut)

	# print the atomic level information
	if parameters.printAtoms:
		atoms.show_info(parameters)

	# print unit cell info, vectors, angles, symmetry ...
	if parameters.printCell:
		lattice.show_info(parameters,cell,parameters.symprec)

	# print bonds and optionally neighbor info
	if parameters.printBonds and parameters.printNlist:
		nn.show_info( parameters,atoms,2 )
	elif parameters.printBonds:
		nn.show_info( parameters,atoms )

	# attempt to reduce convetional cell to primitive cell
	if parameters.getPrimitive:
		primitive = spglib.find_primitive( cell,parameters.symprec )		
		patoms = Atoms(lattice=primitive[0],positions=primitive[1],numbers=primitive[2],fractional=True)
		try:
			test = primitive[0][0][0] * 0
		except TypeError:
			print("Warning: could not reduce to primitive cell...")

		reader.write_vasp( primitive[0],patoms,parameters )


if __name__ == "__main__": main()

