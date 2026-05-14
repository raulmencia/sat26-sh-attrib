# attribution.py

# Authors:
    # Pablo Martinez-Naredo
    # Raul Mencia
    # Joao Marques-Silva
    # Carlos Mencia

from sh_common import *
import gc



class Config:

    def __init__(self, x = None):
        
        self.f_file = None
        self.seed = None
        self.solver_name = "m22"
        self.output_file = None
        
        self.mode = 1
        
        self.original_formula_file = None
        self.MUS_fml_file = None
        self.MCS_fml_file = None
        
        self.iterations = 1000000000
        self.alpha = None
        self.epsilon = None
        
        self.time_limit = 1000000000
        self.timeout_sh = None
        

        if x is not None:
            self.set_config(x)
            self.check_config(x)
    
    
    def set_config(self, x):
        if len(x) == 1:
            print(self.helpString(x))
            exit()
            
        self.mode = int(x[1])
        
        i = 2
        while i < len(x):
            if x[i] == "-seed":
                i += 1
                self.seed = int(x[i])
            elif x[i] == "-sn":
                i += 1
                self.solver_name = x[i]
            elif sys.argv[i] == "-output":
                i += 1
                self.output_file = x[i]

            elif sys.argv[i] == "-it":
                i += 1
                self.iterations = int(x[i])
            elif sys.argv[i] == "-alpha":
                i += 1
                self.alpha = float(x[i])
            elif sys.argv[i] == "-epsilon":
                i += 1
                self.epsilon = float(x[i])
            
            elif sys.argv[i] == "-tl":
                i += 1
                self.time_limit = int(x[i])
            elif x[i] == "-timeout":
                i += 1
                self.timeout_sh = int(x[i])
            
            elif x[i] == "-f":
                i += 1
                self.original_formula_file = x[i]
            elif x[i] == "-mus":
                i += 1
                self.MUS_fml_file = x[i]
            elif x[i] == "-mcs":
                i += 1
                self.MCS_fml_file = x[i]

            else:
                self.error_message("Error: " + str(x[i]) + " is not recognized as a valid argument.", x)

            i += 1


    def check_config(self, x):
        # Checking that the configuration is coherent
        if self.mode == 1:
            if self.original_formula_file is None:
                self.error_message("Error: In mode 1, the original formula (-f) must be provided.", x)
        elif self.mode == 2:
            if self.original_formula_file is None:
                self.error_message("Error: In mode 2, the original formula (-f) must be provided.", x)
            if self.MUS_fml_file is None:
                self.error_message("Error: In mode 2, the MUS formula (-mus) must be provided.", x)
            if self.MCS_fml_file is None:
                self.error_message("Error: In mode 2, the MCS formula (-mcs) must be provided.", x)
        else:
            self.error_message("Error: mode must be 1 or 2.", x)
        
        if self.alpha is not None or self.epsilon is not None:
            if self.alpha is None:
                self.error_message("Error: If epsilon is provided, alpha must be provided as well (-alpha)", x)
            elif self.epsilon is None:
                self.error_message("Error: If alpha is provided, epsilon must be provided as well (-epsilon)", x)
            else:
                self.iterations = get_iterations(self.alpha, self.epsilon)
                print("Target: " + str(self.iterations) + " iterations")
        

    def helpString(self, x):
        s = ""
        s += "Usage: python " + x[0] + " mode arg1 arg2 arg3 ...\n"
        s += "\tmode: used to choose the Shapley-Shubik attribution method (mandatory). Options:\n"
        s += "\t\t1: SAT-based attribution. This requires the following to be provided:\n"
        s += "\t\t\t-f <filename>: filename contains the path to the file with the original formula\n"
        s += "\t\t2: Interval-based attribution. This requires the following to be provided:\n"
        s += "\t\t\t-f <filename>: filename contains the path to the file with the original formula\n"
        s += "\t\t\t-mus <filename>: filename contains the path to the file with the MUS formula\n"
        s += "\t\t\t-mcs <filename>: filename contains the path to the file with the MCS formula\n\n"
        
        s += "\tValid arguments:\n"
        s += "\t\t-seed <n>: n is the seed.\n"
        s += "\t\t-output <filename>: filename is the name of the output file (default: print on console)\n"
        s += "\t\t-sn <solver_name>: solver_name is pysats solver name (default: m22)\n"

        s += "\t\t-it <n>: n is the number of iterations (default: 1000000000)\n"
        s += "\t\tAlternative to -it: Providing alpha and epsilon\n"
        s += "\t\t\t-alpha <alpha>\n"
        s += "\t\t\t-epsilon <epsilon>\n"
        s += "\t\t\tNote: If you use alpha or epsilon, you must use both\n"
        s += "\t\t-tl <s>: s is the time limit in seconds (default: 1000000000)\n"
        s += "\t\t-timeout <s>: s is the TIMEOUT limit in seconds. Requires a MiniSat-like solver (m22, mgh, ...)\n"

        s += "\n\n"
        s += "Examples:\n"
        s += "To run the SAT-based attribution approach (deletion with unsat cores):\n"
        s += "\tpython " + x[0] + " 1 -f ex/ex.cnf -it 960365 -tl 7200 -seed 1\n"
        s += "To run the interval-based attribution approach:\n"
        s += "\tpython " + x[0] + " 2 -f ex/ex.cnf -mus ex/exmus.cnf -mcs ex/exmcs.cnf -it 960365 -tl 7200 -seed 1\n"

        return s
    
    
    def error_message(self, s, x):
        print(s)
        print(self.helpString(x))
        exit()



class ShapleyShubikCalculator:

    def __init__(self, config):
        self.config = config
        self.stats = Stats(config)
        self.tc = TimeControl(self.config.time_limit)
        self.sv = None
        self.interrupted = False


    def compute_indices(self):
        if self.config.mode == 1:
            return self.compute_indices_original()
        else:
            return self.compute_indices_mcs()


    def compute_indices_original(self):

        self.stats.start_time = time.time()
        
        sh_table = {}
        
        fml = OrigFormula(self.config.original_formula_file)
        self.sv = IncrementalSolver(fml, self.config)

        c = CriticalClauseFinder(self.sv, self.config)
        
        order = []

        for i in range(len(self.sv.sel_variables)):
            order.append(i)
            sh_table[self.sv.vC.dv[self.sv.sel_variables[i]]] = 0

        it = 0
        while not self.termination_condition(it):
            random.shuffle(order)
            
            cc = c.find_critical_clause(order)

            if cc is not None:
                sh_table[self.sv.vC.dv[self.sv.sel_variables[order[cc]]]] += 1
                it += 1

        if it > 0:
            for key in sh_table:
                sh_table[key] /= it
        else:
            for key in sh_table:
                sh_table[key] = "-"

        self.sv.delete()

        self.stats.end_time = time.time()
        self.stats.iterations_completed = it
            
        self.output(sh_table)


    def compute_indices_mcs(self):

        self.stats.start_time = time.time()
        
        mus_fml = OrigFormula(self.config.MUS_fml_file)
        self.sv = IncrementalSolver(mus_fml, self.config)
        
        mcs_fml = OrigFormula(self.config.MCS_fml_file)
        sv_mcs = IncrementalSolver(mcs_fml, self.config)
        
        ori_fml = OrigFormula(self.config.original_formula_file)
        ofml_nclauses = len(ori_fml.clauses)
        del ori_fml
        gc.collect()

        c = CriticalClauseFinder(self.sv, self.config)
        
        mcs_ccf = MCS_CCF(sv_mcs, self.config)

        order = []
        for i in range(ofml_nclauses):
            order.append(i)
        
        all_ids = []
        for x in order:
            all_ids.append(x+1)

        osvm = OrderSelVarsManager(self.sv, sv_mcs, order, ofml_nclauses)

        res_mus_mcs = ResultsMUSMCS(self.config, all_ids)

        order_mus = []
        for i in range(len(self.sv.sel_variables)):
            order_mus.append(i)
        
        order_mcs = []
        for i in range(len(sv_mcs.sel_variables)):
            order_mcs.append(i)


        it = 0
        while not self.termination_condition(it):
            random.shuffle(order)
            osvm.prepare_selvars()
            
            cc = c.find_critical_clause(order_mus)

            ccmcs = 0
            if len(mcs_fml.clauses) > 0:
                ccmus = osvm.indexMUStoMCS(cc)
                ccmcs = mcs_ccf.find_critical_element_models(order_mcs, ccmus)
                ccmcs = osvm.indexMCStoORDER(ccmcs)
            
            cc = osvm.indexMUStoORDER(cc)

            if cc is not None and ccmcs is not None:
                res_mus_mcs.update(order, cc, ccmcs)
                it += 1

        self.sv.delete()
        sv_mcs.delete()

        self.stats.end_time = time.time()
        self.stats.iterations_completed = it

        self.output_mus_mcs(res_mus_mcs)


    def termination_condition(self, n_it):
        # Returns True if we have met a termination condition, and False otherwise
        # Termination conditions are:
        #   - Iterations
        #   - Time
        #   - Interruption
        if self.interrupted:
            return True
        if n_it == self.config.iterations:
            return True
        if self.tc.check_time():
            return True
        return False
    
    
    def output(self, shdict = None):
        s = "Shapley-Shubik:\n"
        
        s += "fml: " + self.config.original_formula_file + ";"
        s += " mode: " + str(self.config.mode) + ";"
        s += " seed: " + str(self.config.seed) + ";\n"

        s += self.stats.build_output()
        s += "\n"
        
        if shdict is not None:
            sorted_shdict = dict(sorted(shdict.items(), key=lambda item: item[1], reverse=True))

            s += str(sorted_shdict)
            s += "\n"
        
        if self.config.output_file is not None:
            f = open(self.config.output_file, "a")
            f.write(s)
            f.close()
        else:
            print(s)


    def output_mus_mcs(self, res_mus_mcs = None):
        s = "Shapley-Shubik:\n"
        
        s += "fml: " + self.config.MUS_fml_file + ";"
        s += " mode: " + str(self.config.mode) + ";"
        s += " seed: " + str(self.config.seed) + ";\n"

        s += self.stats.build_output()
        s += "\n"
        
        if res_mus_mcs is not None:
            s += res_mus_mcs.output_string_it(self.stats.iterations_completed)
        
        if self.config.output_file is not None:
            f = open(self.config.output_file, "a")
            f.write(s)
            f.close()
        else:
            print(s)


    def interrupt(self):
        self.sv.solver.interrupt()
        self.interrupted = True



class CriticalClauseFinder:

    def __init__(self, sv, config):
        # self.solver: Incremental solver to be used (with the considered instance)
        self.solver = sv
        self.config = config
        self.dicAsm = {}


    def find_critical_clause(self, order):
        # It returns the index of the critical clause in order
        cc = None
        
        cc = self.Deletion(order)
        
        if cc is None:
            return None
        else:
            return cc


    def Deletion(self, order):
        finished = False
        i = len(order) - 1
        asm = []

        for index in range(len(order)):
            asm.append(self.solver.sel_variables[order[index]])
            self.dicAsm[self.solver.sel_variables[order[index]]] = index

        while not finished:
            while len(asm) - 1 > i:
                asm.pop()

            sat = self.solver.solve(asm)

            if sat is None:
                return None
                
            if sat:
                i += 1

                finished = True
            else:
                i = self.exploit_core()
                i -= 1

        return i


    def exploit_core(self):
        # Returns the index of the last item that appears in the core
        core = self.solver.get_core()
        m = 0
        for x in core:
            if self.dicAsm[x] > m:
                m = self.dicAsm[x]
        return m



class MCS_CCF:

    def __init__(self, sv, config):
        # self.solver: Incremental solver to be used (with the considered instance)
        self.solver = sv
        self.config = config
    
    
    def find_critical_element_models(self, order, cc = None):
        # Returns index in order that is critical
        # based in insertion
        if cc is None:
            return self.deletion(order)

        ccmcs = len(order)
        
        asm = []
        
        while (ccmcs - 1) > cc:
            ccmcs -= 1
            asm.append(self.solver.sel_variables[order[ccmcs]])
        
        
        sat = True
        while sat:
            ccmcs -= 1
            asm.append(self.solver.sel_variables[order[ccmcs]])
            
            sat = self.solver.solve(asm)
            
            if sat:
                ccmcs_new = self.exploit_model(order, ccmcs)
                while ccmcs > ccmcs_new:
                    ccmcs -= 1
                    asm.append(self.solver.sel_variables[order[ccmcs]])

            if sat is None:
                return None
        
        return ccmcs


    def exploit_model(self, order, ccmcs):
        # returns the last index that is sat by the model
        
        model = self.solver.get_model()
        
        finished = False
        while not finished:
            ccmcs -= 1
            
            if ccmcs == -1:
                finished = True
            else:
                x = self.solver.sel_variables[order[ccmcs]]
                
                index = 0
                while index < len(model) and abs(model[index]) != x:
                    index += 1
                
                if index < len(model):
                    if x != model[index]:
                        finished = True
                else:
                    finished = True
            
        ccmcs += 1
        
        return ccmcs



class Sh_interval_updater:
    
    def __init__(self, strategy, factor, sh_table, all_ids):
        
        # self.sh_update_strategy:
            # - 1: Uniform
            # - 2: Exponential
        self.sh_update_strategy = strategy
        self.sh_update_factor = factor

        self.sh_table = sh_table

        self.all_ids = all_ids

        # self.series and self.sumseries are auxiliary arrays
        self.series = [1]
        self.sumseries = [1]


    def update_series(self, objlen):
        while len(self.series) < objlen:
            n = self.series[len(self.series)-1] * self.sh_update_factor
            self.series.append(n)
            self.sumseries.append(self.sumseries[len(self.sumseries)-1]+n)


    def update(self, order, cc, ccmss):
        if self.sh_update_strategy == 1:
            self.uniform(order, cc, ccmss)
        elif self.sh_update_strategy == 2:
            self.distribute(order, cc, ccmss)


    def uniform(self, order, cc, ccmcs):
        i = ccmcs
        while i <= cc:
            self.sh_table[self.all_ids[order[i]]] += (1/(cc-ccmcs+1))
            i += 1


    def distribute(self, order, cc, ccmcs):
        d = cc - ccmcs
        if len(self.series) < (d + 1):
            self.update_series(d+1)
        
        for i in range(d+1):
            self.sh_table[self.all_ids[order[ccmcs+i]]] += (self.series[i]/self.sumseries[d])



class OrderSelVarsManager:
    def __init__(self, sv_mus, sv_mcs, order, ofml_nclauses):
        self.sv_mus = sv_mus
        self.sv_mcs = sv_mcs
        self.ofml_nclauses = ofml_nclauses
        
        # self.order[i] represents clause i + 1
        self.order = order
        
        # self.selvarMUSfml[i] is -1 if clause i + 1 does not appear in the MUS formula;
        #                      If it does, it will be equal to the selector variable in it.
        self.selvarMUSfml = []
        
        # self.selvarMCSfml[i] is -1 if clause i + 1 does not appear in the MCS formula;
        #                      If it does, it will be equal to the selector variable in it.
        self.selvarMCSfml = []
        
        # self.MUStoORDER[i] is de position in self.order of self.sv_mus.sel_varibles[i]
        self.MUStoORDER = []
        
        # self.MCStoORDER[i] is de position in self.order of self.sv_mcs.sel_varibles[i]
        self.MCStoORDER = []
        
        # self.MUStoMCS[i] contains the position of selvar i in MCS (or -1 if it is not present)
        self.MUStoMCS = []
        
        self.initialize()


    def initialize(self):
        for i in range(self.ofml_nclauses):
            # self.order.append(i)
            self.selvarMUSfml.append(-1)
            self.selvarMCSfml.append(-1)
        
        for x in self.sv_mus.sel_variables:
            self.selvarMUSfml[self.sv_mus.vC.dv[x] - 1] = x
            self.MUStoORDER.append(-1)
            self.MUStoMCS.append(-1)
        
        for x in self.sv_mcs.sel_variables:
            self.selvarMCSfml[self.sv_mcs.vC.dv[x] - 1] = x
            self.MCStoORDER.append(-1)


    def prepare_selvars(self):
        # self.order has just been reshuffled
        # self.sv_mus.sel_varibles have to follow the same relative order
        # self.MUStoORDER has to be updated as well
        index_mus = 0
        index_mcs = 0
        
        i = 0
        while i < len(self.order) and index_mus < len(self.sv_mus.sel_variables):

            if self.selvarMUSfml[self.order[i]] != -1:
                self.sv_mus.sel_variables[index_mus] = self.selvarMUSfml[self.order[i]]
                self.MUStoORDER[index_mus] = i
                
                # Note: This works because all the variables in the MCS formula 
                # appear in the MUS formula.
                if self.selvarMCSfml[self.order[i]] != -1:
                    self.sv_mcs.sel_variables[index_mcs] = self.selvarMCSfml[self.order[i]]
                    self.MCStoORDER[index_mcs] = i
                    
                    self.MUStoMCS[index_mus] = index_mcs
                    
                    index_mus += 1
                    index_mcs += 1
                else:
                    self.MUStoMCS[index_mus] = -1
                    index_mus += 1
            i += 1
        
        

    def indexMUStoMCS(self, index):
        if self.MUStoMCS[index] != -1:
            return self.MUStoMCS[index]
        else:
            # Look to the first one to the left (if there is none, return 0)
            x = 0
            while index >= 0 and self.MUStoMCS[index] == -1:
                index -= 1
            if index >= 0:
                x = self.MUStoMCS[index]
            return x


    def indexMUStoORDER(self, index):
        return self.MUStoORDER[index]


    def indexMCStoORDER(self, index):
        return self.MCStoORDER[index]



class ResultsMUSMCS:
    def __init__(self, config, all_ids):
        self.config = config
        self.all_ids = all_ids

        self.low = {}
        self.up = {}
        self.conservative = {}
        self.uniform = {}
        self.exponential = {}

        for x in self.all_ids:
            self.low[x] = 0
            self.up[x] = 0
            self.conservative[x] = 0
            self.uniform[x] = 0
            self.exponential[x] = 0

        self.iupdater_uniform = Sh_interval_updater(1, 0, self.uniform, self.all_ids)
        self.iupdater_exponential = Sh_interval_updater(2, 2, self.exponential, self.all_ids)


    def update(self, order, ccmus, ccmcs):
        # low
        if ccmcs == ccmus:
            self.low[self.all_ids[order[ccmus]]] += 1

        # up
        for i in range(ccmcs, ccmus+1):
            self.up[self.all_ids[order[i]]] += 1
            
        # uniform
        self.iupdater_uniform.update(order, ccmus, ccmcs)
            
        # exponential
        self.iupdater_exponential.update(order, ccmus, ccmcs)
        
        # conservative
        self.conservative[self.all_ids[order[ccmus]]] += 1


    def output_string_it(self, it):
        s = ""

        aux = {}
        for key in self.low:
            aux[key] = self.low[key]/it
        s += "low: " + str(aux) + "\n"
        
        aux.clear()
        for key in self.up:
            aux[key] = self.up[key]/it
        s += "up: " + str(aux) + "\n"
        
        aux.clear()
        for key in self.conservative:
            aux[key] = self.conservative[key]/it
        s += "conservative: " + str(aux) + "\n"
        
        aux.clear()
        for key in self.uniform:
            aux[key] = self.uniform[key]/it
        s += "uniform: " + str(aux) + "\n"
        
        aux.clear()
        for key in self.exponential:
            aux[key] = self.exponential[key]/it
        s += "exponential: " + str(aux) + "\n"

        return s



class Stats:

    def __init__(self, config):
        self.config = config
        self.start_time = None
        self.end_time = None
        self.iterations_completed = 0


    def build_output(self):
        s = ""
        s += "iterations completed: " + str(self.iterations_completed) + "; "
        s += "time: " + str(self.end_time - self.start_time) + ";"
        if self.config.alpha is not None:
            s += " error (alpha = " + str(self.config.alpha) + "): " + str(get_error(self.config.alpha, self.iterations_completed)) + ";"
        return s




config = Config(sys.argv)

if config.seed is not None:
    random.seed(config.seed)

ssc = ShapleyShubikCalculator(config)

timer = None
if config.timeout_sh is not None:
    timer = Timer(config.timeout_sh, interrupt, [ssc])
    timer.start()

ssc.compute_indices()
        
if timer is not None:
    timer.cancel()
