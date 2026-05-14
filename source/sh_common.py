# sh_common.py

# Authors:
    # Pablo Martinez-Naredo
    # Raul Mencia
    # Joao Marques-Silva
    # Carlos Mencia

import random
import time
import sys
import math
import statistics
import lzma

from pysat.solvers import Solver, Minisat22

from threading import Timer



class OrigFormula:

    def __init__(self, fn = None):
        # self.clauses: the soft clauses
        self.clauses = []
        # self.hard_clauses: the hard clauses
        self.hard_clauses = []
        if fn is not None:
            self.formula_from_file(fn)


    def add_clause(self, clause):
        # To add a soft clause
        self.clauses.append(clause)


    def check_cnf_wcnf(self, filename):
        # Returns True if filename contains a formula in CNF format
            # or False if the formula is in WCNF format
        f = get_readfile_object(filename)
        a = f.readline().split()
        f.close()
        
        if a[1] == "Standarized":
            # wcnf
            return False
        else:
            # cnf
            return True


    def formula_from_file(self, fname):
        
        if self.check_cnf_wcnf(fname):
            self.cnf_formula_from_file(fname)
        else:
            self.wcnf_formula_from_file(fname)


    def cnf_formula_from_file(self, fname):
        f = get_readfile_object(fname)
        
        a = f.readline()
        
        while a[0] != "p":
            a = f.readline()

        nclauses = int(a.split(" ")[3])
        
        i = 0
        while i < nclauses:
            x = f.readline().split()
            cs = []
            j = 0
            while x[j][0] != "0":
                cs.append(int(x[j]))
                j += 1
            self.clauses.append(cs)
            i += 1
            
        f.close()
    
    
    def wcnf_formula_from_file(self, fname):
        f = get_readfile_object(fname)
        
        a = f.readline()
        
        while a[0] == "c":
            a = f.readline()
        
        x = a.split()
        
        while len(x) > 1:
            cs = []
            index = 1
            while x[index][0] != "0":
                cs.append(int(x[index]))
                index += 1
            
            if x[0][0] == "h":
                # hard
                self.hard_clauses.append(cs)
            else:
                # soft
                self.clauses.append(cs)
            
            a = f.readline()
            x = a.split()

        f.close()


    def formula_to_file(self, fname):
        # Creates a file called 'fname' and writes a formula on CNF format on it (used to write MUS/MCS formulas)
        variables = set()
        for x in self.clauses:
            for y in x:
                variables.add(abs(y))
                
        f = open(fname, "w")
                
        f.write("p cnf " + str(len(variables)) + " " + str(len(self.clauses)) + "\n")
        for x in self.clauses:
            for y in x:
                f.write(str(y) + " ")
            f.write("0\n")

        f.close()



class VariablesControl:

    def __init__(self):
        # self.n_variables contains the number of considered variables (formula + selector variables)
        self.n_variables = 0
        # self.variables contains the variables in this structure associated to those in the formula.
        self.variables = {}
        # self.sel contains the selector variables associated with the clauses
        self.sel = {}
        # self.dv is a dictionary associating self.selvm entries with original clauses
        self.dv = {}


    def get_var(self, x):
        # If there is an entry in self.variables for x, it returns it.
        # Else:
        #   Creates and returns a new variable K associated to x
        #       - self.variables[x] = K
        #   Adds a new entry to the self.dv
        #       - self.dv[K] = x
        v = abs(x)
        if v not in self.variables:
            self.n_variables += 1
            self.variables[v] = self.n_variables
            self.dv[self.n_variables] = v
        return self.variables[v]


    def get_sel(self, x):
        # If there is an entry in self.sel for x, it returns it.
        # Else:
        #   Creates and returns a new variable K associated to x
        #       - self.sel[x] = K
        #   Adds a new entry to the self.dv
        #       - self.dv[K] = x
        v = abs(x)
        if v not in self.sel:
            self.n_variables += 1
            self.sel[v] = self.n_variables
            self.dv[self.n_variables] = v
        return self.sel[v]
    
    
    def get_original(self, x):
        # Returns the original value (selvar / variable) associated with x
        return self.dv[x]



class IncrementalSolver:

    def __init__(self, orig_formula, config):
        # self.orig_formula: Object of class 'OrigFormula' containing the original formula
        self.orig_formula = orig_formula
        # self.config: Configuration object of class 'Config'
        self.config = config
        
        #self.solver: Pysat solver (self.config.solver_name contains the pysats solver name)
        self.solver = Solver(name=self.config.solver_name)
        
        # self.vC: Object of class 'VariablesControl', used to manage variables
        self.vC = VariablesControl()
        # self.sel_variables: List with the selector variables
        self.sel_variables = []
        
        self.initialize()


    def initialize(self):
        # Initializes the incremental solvers structures:
        #   - self.solver: clauses are added to it
        #   - self.vC
        #   - self.sel_variables
        
        if self.config.mode == 1:
            # In 'clause mode', each clause has its selector variable
            for i in range(len(self.orig_formula.clauses)):
                # soft
                x = self.orig_formula.clauses[i]
                
                solver_clause = []
                
                sv = self.vC.get_sel(i+1)
                self.sel_variables.append(sv)
                solver_clause.append(sv*-1)
                
                for v in x:
                    if v < 0:
                        solver_clause.append(self.vC.get_var(v)*-1)
                    else:
                        solver_clause.append(self.vC.get_var(v))

                self.solver.add_clause(solver_clause)
            
            for i in range(len(self.orig_formula.hard_clauses)):
                # hard
                x = self.orig_formula.hard_clauses[i]
                
                solver_clause = []
                
                for v in x:
                    if v < 0:
                        solver_clause.append(self.vC.get_var(v)*-1)
                    else:
                        solver_clause.append(self.vC.get_var(v))

                self.solver.add_clause(solver_clause)
                
        elif self.config.mode == 2:
            # In MUS mode, we don't create additional selector variables
            # All the variables are selector variables
            for x in self.orig_formula.clauses:
                solver_clause = []
                
                for v in x:
                    solver_clause.append(self.vC.get_sel(v)*-1)

                self.solver.add_clause(solver_clause)
                
            for slv in self.vC.sel:
                self.sel_variables.append(self.vC.sel[slv])


    def solve(self, asmps = None, phs = None):
        if asmps is None:
            asmps = self.sel_variables
        if phs is not None:
            self.solver.set_phases(phs)
        return self.solver.solve_limited(assumptions=asmps, expect_interrupt=True)


    def get_model(self):
        return self.solver.get_model()


    def get_core(self):
        return self.solver.get_core()


    def delete(self):
        self.solver.delete()



class TimeControl:

    def __init__(self, tl):
        self.time_limit = tl
        self.start_time = time.time()


    def check_time(self):
        # Returns True if time has run out
        return ((time.time() - self.start_time) > self.time_limit)



def get_iterations(alpha, epsilon):
    alpha_half = alpha / 2
    sigma_square = 0.25
    dist = statistics.NormalDist(mu=0, sigma=1)
    X = dist.inv_cdf(alpha_half)
    m = X * X * sigma_square / (epsilon * epsilon)
    
    return int(math.ceil(m))



def get_error(alpha, iterations):
    alpha_half = alpha / 2
    sigma_square = 0.25
    dist = statistics.NormalDist(mu=0, sigma=1)
    X = dist.inv_cdf(alpha_half)
    error = math.sqrt(X * X * sigma_square / iterations)

    return error



def is_xz(filename):
    try:
        with open(filename, "rb") as f:
            s = f.read(6)
            return s == b'\xfd7zXZ\x00'
    except IOError:
        return False



def get_readfile_object(filename):
    # Returns a file object associated with file 'filename'.
    # This is used to check if the file is in xz first.
    f = None
        
    if is_xz(filename):
        f = lzma.open(filename, "rt")
    else:
        f = open(filename, "r")

    return f



def interrupt(x):
    x.interrupt();

