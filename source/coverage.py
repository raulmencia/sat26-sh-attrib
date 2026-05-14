# coverage.py

# Authors:
    # Pablo Martinez-Naredo
    # Raul Mencia
    # Joao Marques-Silva
    # Carlos Mencia

from pysat.allies.approxmc import Counter
from pysat.formula import CNF
import sys
from decimal import Decimal

class CoverageCalculator:

    def __init__(self, config):

        self.config = config
        
        self.variables = []
        
        self.cnf = CNF()
        
        mus_cs = self.get_clauses_from_file(self.config.filename_mus)
        self.append_mus_clauses(mus_cs)
        
        if self.config.filename_mcs is not None:
            mcs_cs = self.get_clauses_from_file(self.config.filename_mcs)
            self.append_mcs_clauses(mcs_cs)


    def n_models(self):
        c = Counter(self.cnf, epsilon=self.config.epsilon, delta=self.config.delta)
        return c.count()


    def get_coverage(self):
        topv = max(self.variables)
        nv = topv + 1
        self.cnf.append([nv])

        x = Decimal(self.n_models())
        n_m = x
                
        for i in range(topv):
            x /= 2
        
        cv = 1 - x

        self.output(str(cv))


    def output(self, s):
        if self.config.output_file is not None:
            f = open(self.config.output_file, "w")
            f.write(s)
            f.close()
        else:
            print(s)


    def append_mus_clauses(self, clauses):
        for c in clauses:
            for x in c:
                v = abs(x)
                if v not in self.variables:
                    self.variables.append(v)
            self.cnf.append(c)


    def append_mcs_clauses(self, clauses):
        for c in clauses:
            nc = []
            for x in c:
                v = abs(x)
                nc.append(v)
                if v not in self.variables:
                    self.variables.append(v)
            self.cnf.append(nc)


    def get_clauses_from_file(self, filename):
        f = open(filename, 'r')
        
        f.readline()

        ncs = []
        
        for line in f:
            x = line.split()
            
            if len(x) > 0:
                cls = []
                
                index = 0
                
                while x[index][0] != "0":
                    cls.append(int(x[index]))
                    index += 1
                
                ncs.append(cls)
                
        f.close()
        return ncs



class Config:

    def __init__(self, x = None):
        self.filename_mus = ""
        self.filename_mcs = None
        self.output_file = None
        self.delta = 0.95
        self.epsilon = 0.025
        
        if x is not None:
            self.set_config(x)


    def set_config(self, x):
        if len(x) == 1:
            print(self.helpString(x))
            exit()
            
        self.filename_mus = x[1]
        
        i = 2
        while i < len(x):
            if sys.argv[i] == "-mcs":
                i += 1
                self.filename_mcs = x[i]
                
            elif sys.argv[i] == "-output":
                i += 1
                self.output_file = x[i]

            elif sys.argv[i] == "-d":
                i += 1
                self.delta = float(x[i])

            elif sys.argv[i] == "-e":
                i += 1
                self.epsilon = float(x[i])

            else:
                print("Error: " + str(x[i]) + " is not recognized as a valid argument.")
                print(self.helpString(x))
                exit()
            i += 1


    def helpString(self, x):
        s = ""
        
        s += "Usage: python " + x[0] + " filename_mus arg1 arg2 arg3 ...\n"
        s += "\tfilename_mus: contains the path to the file with the MUS formula. (Mandatory)\n"
        s += "\t-output <filename>: filename is the name of the output file (default: print on console)\n"
        s += "\t-mcs <filename_mcs>: filename_mcs contains the path to the file with the MUS formula.\n"
        s += "\t-e <epsilon> (default: 0.025)\n"
        s += "\t-d <delta> (default: 0.95)\n"
        
        return s


config = Config(sys.argv)

cov = CoverageCalculator(config)

cov.get_coverage()






