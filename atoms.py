# -*- coding: utf-8 -*-

import numpy as np
import sys

class Atoms:
    def __init__(self,
                 lattice=None,
                 volume=None,
                 symbols=None,
                 positions=None,
                 numbers=None, 
                 ntypes=None,
                 fractional=True):

        # positions (cartesian)
        self.x = None
        # scaled positions (lattice)
        self.xs = None

        self.unq = None

        self.v = volume

        if not ntypes is None:
            self.ntypes = ntypes
        else:
            self.ntypes = []

        self.ids = []
        if not symbols is None:
            self.set_ids(symbols)
    
        # matrix for cell vectors
        self.H = []
        for i in range(3):
            ax = []
            ax.append(lattice[i][0])
            ax.append(lattice[i][1])
            ax.append(lattice[i][2])
            self.H.append( ax )

        if fractional:
            self.xs = positions
            self.x = self.cart_coordinates(lattice)
        else:
            self.x = positions
            self.xs = self.reduce_coordinates(lattice)

        # Atom symbols
        self.symbols = symbols

        # Atomic numbers
        if numbers is None:
            self.numbers = None
        else:
            self.numbers = np.array(numbers, dtype='intc')

        # masses
        self.masses = None
        self.total_mass =0

        # number --> symbol
        if not self.numbers is None:
            self._numbers_to_symbols()

        # symbol --> number
        elif not self.symbols is None:
            self._symbols_to_numbers()

        # symbol --> mass
        if self.symbols and (self.masses is None):
            self._symbols_to_masses()

        ## finish construction

    def unique_items(self,seq):
        # order preserving
        noDupes = []
        [noDupes.append(i) for i in seq if not noDupes.count(i)]
        return noDupes

    def count_types(self,ids):
        # ntypes = []

        n = self.unique_items(ids)

        for i in range(len(n)):
            self.ntypes.append(0)

        count = 1
        t = 0

        for i in range(0,len(ids)-1):
            tmp = ids[i]
            if tmp == ids[i+1]:
                count += 1
            else:
                self.ntypes[t] = count
                t += 1
                count = 1
            
        if ids[len(ids)-2] != ids[len(ids)-1]:
            self.ntypes[t] = 1
        else:
            self.ntypes[t] = count

        # print "NT", len(self.ntypes)
        # print self.ntypes

        return self.ntypes


    def get_compound( self ):
        self.unq = self.unique_items(self.symbols)
        name = self.unq[0]
        name += str(self.ntypes[0])

        for i in range(1,len(self.unq)):
            name += self.unq[i]
            name += str(self.ntypes[i])

        return name

    def set_ids(self,symbols):        
        n = 0
        # print len(symbols)
        for i in range(0,len(symbols)-1):
            n += 1
            atomi = symbols[i]
            atomj = symbols[i+1]
            if atomi == atomj:
                self.ids.append(n)
            else:
                self.ids.append(n)
                n = 0
        if symbols[len(symbols)-2] != symbols[len(symbols)-1] :
            self.ids.append(1)
        else:
            self.ids.append(n)

    def show_info(self,parameters):
        if parameters.save:
            name = self.get_compound() + ".atoms"
            out = open(name,'w')
        else:
            out = sys.stdout

        try:
            if out is sys.stdout:
                out.write('\n' + "/*-- Atoms --*/" + '\n')
            else:
                out.write("/*-- Atoms --*/" + '\n')
            out.write("Structure     = %s    " % parameters.FILE + '\n')
            out.write("Compound      = %s    " % self.get_compound() + '\n')
            out.write("Natoms        = %i    " % self.get_number_of_atoms() + '\n')
            out.write("Ntypes        = %i    " % len(self.ntypes) + '\n')
            out.write("Species       =")
            for i in range(len(self.unq)):
                out.write(" %s" % self.unq[i])
            out.write('\n')
            out.write("Type count    =")
            for i in range(len(self.ntypes)):
                out.write(" %i" % self.ntypes[i])
            out.write('\n')
            out.write('Masses        = ')
            for i in range(len(self.masses)):
                out.write("%f " % self.masses[i])
            out.write('amu\n')
            out.write("Density       = %f g/cm^3\n" % self.get_density())
        finally:
            if parameters.save:
                out.close()

    def show_positions(self):
        print("Lattice coordinates:")
        for i in range(len(self.x)):
            print("%s%i %5f %5f %5f" % (self.symbols[i],self.ids[i],self.xs[i][0],self.xs[i][1],self.xs[i][2]))

        print("Cartesian coordinates:")
        for i in range(len(self.x)):
            print("%s%i %5f %5f %5f" % (self.symbols[i],self.ids[i],self.x[i][0],self.x[i][1],self.x[i][2]))

    def get_pos(self):
        return self.x

    def get_spos(self):
        return self.xs

    def reduce_coordinates(self,lattice):
        return np.dot(self.x, np.linalg.inv(lattice))

    def cart_coordinates(self,lattice):
        return np.dot(self.xs, lattice)

    def set_masses(self, masses):
        if masses is None:
            self.masses = None
        else:
            self.masses = np.array(masses, dtype='double')

    def get_masses(self):
        if self.masses is None:
            return None
        else:
            return self.masses.copy()

    def get_total_mass(self):
        self.total_mass =0
        for i in range(len(self.masses)):
            x = self.ntypes[i]
            m = self.masses[i]
            self.total_mass += x*m
        return self.total_mass

    def get_density(self):
        self.get_total_mass()
        amu2grams = 1.66054e-24
        gmass = self.total_mass * amu2grams
        # totalMass = natoms * gmass
        ang2cc = 1e-24
        volcc = self.v * ang2cc
        return gmass / volcc

    def get_chemical_symbols(self):
        return self.symbols[:]

    def get_number_of_atoms(self):
        return len(self.x)

    def get_atomic_numbers(self):
        return self.numbers.copy()

    def get_volume(self):
        return np.linalg.det(self.cell)
        
    def _numbers_to_symbols(self):
        self.symbols = [atom_data[n][1] for n in self.numbers]
        
    def _symbols_to_numbers(self):
        self.numbers = np.array([symbol_map[s]
                                 for s in self.symbols], dtype='intc')
        
    def _symbols_to_masses(self):
        u = self.unique_items(self.symbols)
        masses = [atom_data[symbol_map[s]][3] for s in u]
        if None in masses:
            self.masses = None
        else:
            self.masses = np.array(masses, dtype='double')

atom_data = [ 
    [  0, "X", "X", None], # 0
    [  1, "H", "Hydrogen", 1.00794], # 1
    [  2, "He", "Helium", 4.002602], # 2
    [  3, "Li", "Lithium", 6.941], # 3
    [  4, "Be", "Beryllium", 9.012182], # 4
    [  5, "B", "Boron", 10.811], # 5
    [  6, "C", "Carbon", 12.0107], # 6
    [  7, "N", "Nitrogen", 14.0067], # 7
    [  8, "O", "Oxygen", 15.9994], # 8
    [  9, "F", "Fluorine", 18.9984032], # 9
    [ 10, "Ne", "Neon", 20.1797], # 10
    [ 11, "Na", "Sodium", 22.98976928], # 11
    [ 12, "Mg", "Magnesium", 24.3050], # 12
    [ 13, "Al", "Aluminium", 26.9815386], # 13
    [ 14, "Si", "Silicon", 28.0855], # 14
    [ 15, "P", "Phosphorus", 30.973762], # 15
    [ 16, "S", "Sulfur", 32.065], # 16
    [ 17, "Cl", "Chlorine", 35.453], # 17
    [ 18, "Ar", "Argon", 39.948], # 18
    [ 19, "K", "Potassium", 39.0983], # 19
    [ 20, "Ca", "Calcium", 40.078], # 20
    [ 21, "Sc", "Scandium", 44.955912], # 21
    [ 22, "Ti", "Titanium", 47.867], # 22
    [ 23, "V", "Vanadium", 50.9415], # 23
    [ 24, "Cr", "Chromium", 51.9961], # 24
    [ 25, "Mn", "Manganese", 54.938045], # 25
    [ 26, "Fe", "Iron", 55.845], # 26
    [ 27, "Co", "Cobalt", 58.933195], # 27
    [ 28, "Ni", "Nickel", 58.6934], # 28
    [ 29, "Cu", "Copper", 63.546], # 29
    [ 30, "Zn", "Zinc", 65.38], # 30
    [ 31, "Ga", "Gallium", 69.723], # 31
    [ 32, "Ge", "Germanium", 72.64], # 32
    [ 33, "As", "Arsenic", 74.92160], # 33
    [ 34, "Se", "Selenium", 78.96], # 34
    [ 35, "Br", "Bromine", 79.904], # 35
    [ 36, "Kr", "Krypton", 83.798], # 36
    [ 37, "Rb", "Rubidium", 85.4678], # 37
    [ 38, "Sr", "Strontium", 87.62], # 38
    [ 39, "Y", "Yttrium", 88.90585], # 39
    [ 40, "Zr", "Zirconium", 91.224], # 40
    [ 41, "Nb", "Niobium", 92.90638], # 41
    [ 42, "Mo", "Molybdenum", 95.96], # 42
    [ 43, "Tc", "Technetium", None], # 43
    [ 44, "Ru", "Ruthenium", 101.07], # 44
    [ 45, "Rh", "Rhodium", 102.90550], # 45
    [ 46, "Pd", "Palladium", 106.42], # 46
    [ 47, "Ag", "Silver", 107.8682], # 47
    [ 48, "Cd", "Cadmium", 112.411], # 48
    [ 49, "In", "Indium", 114.818], # 49
    [ 50, "Sn", "Tin", 118.710], # 50
    [ 51, "Sb", "Antimony", 121.760], # 51
    [ 52, "Te", "Tellurium", 127.60], # 52
    [ 53, "I", "Iodine", 126.90447], # 53
    [ 54, "Xe", "Xenon", 131.293], # 54
    [ 55, "Cs", "Caesium", 132.9054519], # 55
    [ 56, "Ba", "Barium", 137.327], # 56
    [ 57, "La", "Lanthanum", 138.90547], # 57
    [ 58, "Ce", "Cerium", 140.116], # 58
    [ 59, "Pr", "Praseodymium", 140.90765], # 59
    [ 60, "Nd", "Neodymium", 144.242], # 60
    [ 61, "Pm", "Promethium", None], # 61
    [ 62, "Sm", "Samarium", 150.36], # 62
    [ 63, "Eu", "Europium", 151.964], # 63
    [ 64, "Gd", "Gadolinium", 157.25], # 64
    [ 65, "Tb", "Terbium", 158.92535], # 65
    [ 66, "Dy", "Dysprosium", 162.500], # 66
    [ 67, "Ho", "Holmium", 164.93032], # 67
    [ 68, "Er", "Erbium", 167.259], # 68
    [ 69, "Tm", "Thulium", 168.93421], # 69
    [ 70, "Yb", "Ytterbium", 173.054], # 70
    [ 71, "Lu", "Lutetium", 174.9668], # 71
    [ 72, "Hf", "Hafnium", 178.49], # 72
    [ 73, "Ta", "Tantalum", 180.94788], # 73
    [ 74, "W", "Tungsten", 183.84], # 74
    [ 75, "Re", "Rhenium", 186.207], # 75
    [ 76, "Os", "Osmium", 190.23], # 76
    [ 77, "Ir", "Iridium", 192.217], # 77
    [ 78, "Pt", "Platinum", 195.084], # 78
    [ 79, "Au", "Gold", 196.966569], # 79
    [ 80, "Hg", "Mercury", 200.59], # 80
    [ 81, "Tl", "Thallium", 204.3833], # 81
    [ 82, "Pb", "Lead", 207.2], # 82
    [ 83, "Bi", "Bismuth", 208.98040], # 83
    [ 84, "Po", "Polonium", None], # 84
    [ 85, "At", "Astatine", None], # 85
    [ 86, "Rn", "Radon", None], # 86
    [ 87, "Fr", "Francium", None], # 87
    [ 88, "Ra", "Radium", None], # 88
    [ 89, "Ac", "Actinium", None], # 89
    [ 90, "Th", "Thorium", 232.03806], # 90
    [ 91, "Pa", "Protactinium", 231.03588], # 91
    [ 92, "U", "Uranium", 238.02891], # 92
    [ 93, "Np", "Neptunium", None], # 93
    [ 94, "Pu", "Plutonium", None], # 94
    [ 95, "Am", "Americium", None], # 95
    [ 96, "Cm", "Curium", None], # 96
    [ 97, "Bk", "Berkelium", None], # 97
    [ 98, "Cf", "Californium", None], # 98
    [ 99, "Es", "Einsteinium", None], # 99
    [100, "Fm", "Fermium", None], # 100
    [101, "Md", "Mendelevium", None], # 101
    [102, "No", "Nobelium", None], # 102
    [103, "Lr", "Lawrencium", None], # 103
    [104, "Rf", "Rutherfordium", None], # 104
    [105, "Db", "Dubnium", None], # 105
    [106, "Sg", "Seaborgium", None], # 106
    [107, "Bh", "Bohrium", None], # 107
    [108, "Hs", "Hassium", None], # 108
    [109, "Mt", "Meitnerium", None], # 109
    [110, "Ds", "Darmstadtium", None], # 110
    [111, "Rg", "Roentgenium", None], # 111
    [112, "Cn", "Copernicium", None], # 112
    [113, "Uut", "Ununtrium", None], # 113
    [114, "Uuq", "Ununquadium", None], # 114
    [115, "Uup", "Ununpentium", None], # 115
    [116, "Uuh", "Ununhexium", None], # 116
    [117, "Uus", "Ununseptium", None], # 117
    [118, "Uuo", "Ununoctium", None], # 118
    ]

symbol_map = {
    "H":1,
    "He":2,
    "Li":3,
    "Be":4,
    "B":5,
    "C":6,
    "N":7,
    "O":8,
    "F":9,
    "Ne":10,
    "Na":11,
    "Mg":12,
    "Al":13,
    "Si":14,
    "P":15,
    "S":16,
    "Cl":17,
    "Ar":18,
    "K":19,
    "Ca":20,
    "Sc":21,
    "Ti":22,
    "V":23,
    "Cr":24,
    "Mn":25,
    "Fe":26,
    "Co":27,
    "Ni":28,
    "Cu":29,
    "Zn":30,
    "Ga":31,
    "Ge":32,
    "As":33,
    "Se":34,
    "Br":35,
    "Kr":36,
    "Rb":37,
    "Sr":38,
    "Y":39,
    "Zr":40,
    "Nb":41,
    "Mo":42,
    "Tc":43,
    "Ru":44,
    "Rh":45,
    "Pd":46,
    "Ag":47,
    "Cd":48,
    "In":49,
    "Sn":50,
    "Sb":51,
    "Te":52,
    "I":53,
    "Xe":54,
    "Cs":55,
    "Ba":56,
    "La":57,
    "Ce":58,
    "Pr":59,
    "Nd":60,
    "Pm":61,
    "Sm":62,
    "Eu":63,
    "Gd":64,
    "Tb":65,
    "Dy":66,
    "Ho":67,
    "Er":68,
    "Tm":69,
    "Yb":70,
    "Lu":71,
    "Hf":72,
    "Ta":73,
    "W":74,
    "Re":75,
    "Os":76,
    "Ir":77,
    "Pt":78,
    "Au":79,
    "Hg":80,
    "Tl":81,
    "Pb":82,
    "Bi":83,
    "Po":84,
    "At":85,
    "Rn":86,
    "Fr":87,
    "Ra":88,
    "Ac":89,
    "Th":90,
    "Pa":91,
    "U":92,
    "Np":93,
    "Pu":94,
    "Am":95,
    "Cm":96,
    "Bk":97,
    "Cf":98,
    "Es":99,
    "Fm":100,
    "Md":101,
    "No":102,
    "Lr":103,
    "Rf":104,
    "Db":105,
    "Sg":106,
    "Bh":107,
    "Hs":108,
    "Mt":109,
    "Ds":110,
    "Rg":111,
    "Cn":112,
    "Uut":113,
    "Uuq":114,
    "Uup":115,
    "Uuh":116,
    "Uus":117,
    "Uuo":118,
    }
