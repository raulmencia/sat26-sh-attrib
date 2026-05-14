# enumeration.py

# Authors:
    # Pablo Martinez-Naredo
    # Raul Mencia
    # Joao Marques-Silva
    # Carlos Mencia

from sh_common import *

class ConfigMarco:

    def __init__(self, x = None):
        self.mode = 1
        self.finput = ""
        self.solver_name = "m22"
        self.output_file = None

        self.cmusformula = False
        self.foutput = None
        self.cmcsformula = False
        self.mcsf = None
        
        self.iterations = 1000000000
        self.time_limit = 1000000000
        self.nmus = 1000000000

        self.timeout_marco = None
        
        if x is not None:
            self.set_config(x)


    def set_config(self, x):
        if len(x) == 1:
            print(self.helpString(x))
            exit()
            
        self.finput = x[1]
        
        i = 2
        while i < len(x):
            if x[i] == "-it":
                i += 1
                self.iterations = int(x[i])
            elif x[i] == "-tl":
                i += 1
                self.time_limit = int(x[i])
            elif x[i] == "-nmus":
                i += 1
                self.nmus = int(x[i])
                
            elif x[i] == "-sn":
                i += 1
                self.solver_name = x[i]

            elif x[i] == "-ftf":
                self.cmusformula = True
                i += 1
                self.foutput = x[i]
            elif x[i] == "-mcsf":
                self.cmcsformula = True
                i += 1
                self.mcsf = x[i]
            
            elif sys.argv[i] == "-output":
                i += 1
                self.output_file = x[i]

            elif x[i] == "-timeout":
                i += 1
                self.timeout_marco = int(x[i])
            else:
                print("Error: " + str(x[i]) + " is not recognized as a valid argument.")
                print(self.helpString(x))
                exit()
            i += 1


    def helpString(self, x):
        s = ""
        s += "Usage: python " + x[0] + " filename arg1 arg2 arg3 ...\n"
        
        s += "\tfilename: contains the path to the file with the formula. (Mandatory)\n"
        s += "\tValid arguments:\n"
        
        s += "\t-ftf <f_output>: used to write the MUS formula in file f_output\n"
        s += "\t-mcsf <f_output>: used to write the MCS formula in file f_output\n"
        
        s += "\t-output <filename>: filename is the name of the output file (default: print on console)\n"
        
        s += "\t-sn <solver_name>: solver_name is pysats solver name (default: m22)\n"
        
        s += "\t-it <n>: n is the limit for the number of iterations (default: 1000000000)\n"
        s += "\t-tl <s>: s is the time limit in seconds (default: 1000000000)\n"
        s += "\t-nmus <m>: m is the limit for the number of MUSes found(default: 1000000000)\n"
        s += "\t-timeout <s>: s is the TIMEOUT limit in seconds. Warning: use a MiniSat-like solver (m2, mgh, ...)\n\n"

        return s



class MUSEnumerator:

    def __init__(self, config, solver, sel_vars, stats):
        self.config = config
        self.solver = solver;
        
        self.stats = stats

        self.H = Solver(name=self.config.solver_name)
        self.sv = sel_vars;
        self.sv_H = [];
        self.sv_Hset = set();
        
        self.tc = TimeControl(self.config.time_limit)
        
        self.MUS_formula = None
        if self.config.cmusformula:
            self.MUS_formula = OrigFormula()
        
        self.MCS_formula = None
        if self.config.cmcsformula:
            self.MCS_formula = OrigFormula()
        
        # self.full_enum is True only if a full enumeration of MUSes is achieved
        self.full_enum = False
        
        self.interrupted = False



    def add_clause_to_H(self, C):
        for i in C:
            if abs(i) not in self.sv_Hset:
                self.sv_H = self.sv_H + [abs(i)]
                self.sv_Hset.add(abs(i))
                
        self.H.add_clause(C);


    def compute_max_model(self, model):
        M = [];
        U = set();
        
        # add svs in model
        for s in self.sv_H:
            if model[s-1] > 0:
                M = M + [s]
            else:
                U.add(s);

        # vars outside the map are added to M
        for s in self.sv:
            if s not in self.sv_Hset:
                M = M + [s]

        while len(U) > 0:
            current = U.pop();

            self.H.set_phases(literals=list(U));
            result = self.H.solve(M + [current]);

            if result:
                M = M + [current]
                Hmodel = self.H.get_model();
                to_remove = set()
                for i in U:
                    if Hmodel[i-1] > 0:
                        M = M + [i];
                        to_remove.add(i)
                U = U - to_remove;

        return M;


    def compute_mus(self, M):
        R = set(M);
        MU = set();
        while len(R) > 0:
            current = R.pop();

            res = self.solver.solve(list(MU|R))
            if res is None:
                return None
            if res:
                # current is necessary
                MU.add(current)
            else:
                # reduce R to its core (removing clauses in MU)
                R = set(self.solver.get_core())
                R = R - MU
        return sorted(MU);

		
    def block_mcs(self, C):
        self.stats.nmcs += 1
        self.add_clause_to_H(C);
        
        if self.config.cmcsformula:
            ncl = []
            for m in C:
                x = self.solver.vC.get_original(m)
                ncl.append(-x)
            self.MCS_formula.add_clause(ncl)


    def block_mss(self, M):
        
        Mset = set();
        for i in M:
            Mset.add(i)
		
        C = []
        for i in self.sv:
            if i not in Mset:
                C = C + [i]
        C.sort();
        self.block_mcs(C)


    def block_mus(self, M):
        self.stats.nmus += 1

        self.add_clause_to_H([-m for m in M])
        
        if self.config.cmusformula:
            ncl = []
            for m in M:
                x = self.solver.vC.get_original(m)
                ncl.append(-x)
            self.MUS_formula.add_clause(ncl)


    def compute_next(self):
        result_H = self.H.solve();

        if result_H:
            M = self.compute_max_model(self.H.get_model());
            
            result_solver = self.solver.solve(M);

            if result_solver is None:
                return None;
            if result_solver:
                # M induces an MSS
                self.block_mss(M);
            else:
                MU = self.compute_mus(self.solver.get_core());
                if MU is None:
                    return None;
                self.block_mus(MU);

            return True;

        else:
            return False;


    def enumerate(self):
        res = True
        i = 0
        while not self.termination_condition(res, i):
            res = self.compute_next();
            i += 1
        
        self.stats.iterations_completed = i


    def termination_condition(self, res, n_it):
        # Returns True if we have met a termination condition, and False otherwise
        # Termination conditions are:
        #   - Iterations
        #   - Time
        #   - Not res
        #   - Number of MUSes found
        #   - Interruption
        if self.interrupted:
            return True
        if not res:
            self.full_enum = True
            return True
        if n_it == self.config.iterations:
            return True
        if self.tc.check_time():
            return True
        if self.stats.nmus == self.config.nmus:
            return True
        return False



class MarcoLauncher:
    
    def __init__(self, cfm):
        self.config = cfm
        self.stats = StatsMarco(self.config)
        self.me = None
        self.g = None


    def run_marco(self):
        self.stats.start_time = time.time()

        formula = OrigFormula(self.config.finput)
        self.g = IncrementalSolver(formula, self.config)
        
        self.me = MUSEnumerator(self.config, self.g, self.g.sel_variables, self.stats);
        self.me.enumerate();
        
        self.g.delete()
        self.me.H.delete()
        
        self.stats.end_time = time.time()

        if self.config.foutput is not None:
            self.me.MUS_formula.formula_to_file(self.config.foutput)
        
        if self.config.mcsf is not None:
            self.me.MCS_formula.formula_to_file(self.config.mcsf)
        
        self.output()


    def output(self):
        s = "MARCO: "
        if self.me.full_enum:
            s += "Full enumeration of MUSes achieved"
        s += "\n"

        s += self.stats.build_output()
        s += "\n"

        if self.config.output_file is not None:
            f = open(self.config.output_file, "a")
            f.write(s)
            f.close()
        else:
            print(s)


    def output_line(self, s):
        if self.config.output_file is not None:
            s += "\n"
            f = open(self.config.output_file, "a")
            f.write(s)
            f.close()
        else:
            print(s)


    def interrupt(self):
        self.g.solver.interrupt();
        self.me.interrupted = True



class StatsMarco:

    def __init__(self, config):
        self.config = config

        self.start_time = None
        self.end_time = None
        self.nmus = 0;
        self.nmcs = 0;
        self.iterations_completed = 0
    
    
    def build_output(self):
        s = ""
        s += "nmus: " + str(self.nmus) + " "
        s += "nmcs: " + str(self.nmcs) + " "
        s += "iterations completed: " + str(self.iterations_completed) + " "
        s += "time: " + str(self.end_time - self.start_time)
        
        return s



config = ConfigMarco(sys.argv)

ml = MarcoLauncher(config)

timer = None
if config.timeout_marco is not None:
    timer = Timer(config.timeout_marco, interrupt, [ml])
    timer.start()

ml.run_marco()
        
if timer is not None:
    timer.cancel()

