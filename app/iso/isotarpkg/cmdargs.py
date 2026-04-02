# -*- coding: utf-8 -*-

''' Libraries '''
import re
import argparse
import os, errno
import os.path
import sys
import psutil

class CommandLineArgsParser():
        
    def __init__(self):
        
        self.consensus = None
        self.prediction = None
        self.enrichment = None
        self.no_cpu_cores = None
        self.min_consensus = None
        self.pvalue = None
        self.pvalue_type = None
        self.tools = []
        
        # Define the AgumentParser object
        self.cmdLineArgs = argparse.ArgumentParser(description = 
                                                       "Run run.py")
        
        tools_id_helper = '''Specify a list with tools ID. Available tools ID 
are reported as follow: (0) miRanda, (1) miRmap, (2) TargetScan, (3) 
RNAhybrid, (4) PITA. Example: "0,1,4"'''
        
        cpu_cores_helper = "Specify the right number of CPU cores to be used"
        
        consensus_helper = '''Specify how to combine tools output according to 
the following format: '[AND=0,1,4][OR=3], in which 'AND' and 'OR' stand for 
intersection and union, respectively. Such a format intersects the output of 
tools 0, 1, and 4. The resulting intersection is then joined with the one of 
tool 3'''
        
        prediction_helper = '''Specify 'yes' or 'no' for miRNA(s) prediction 
analysis'''
        
        enrichment_helper = '''Specify 'yes' or 'no' for miRNA(s) enrichment 
analysis'''
        
        enrichment_type_helper = '''Specify which type of enrichment 
analysis, as list, you would like to perform [go, kegg, reactome]'''        
        
        min_consensus_helper = '''Specify a minimum consensus [1-5]'''        
                
        pvalue_helper = '''Specify a pvalue'''
        
        pvalue_type_helper = '''Specify the pvalue type [raw, adj]'''
        
        # No. of cores command line argument:
        self.cmdLineArgs.add_argument('--cores', 
                                      type = int, 
                                      default = 1, 
                                      help = cpu_cores_helper) 
        
        # Tools list command line argument:
        self.cmdLineArgs.add_argument('--tools', 
                                      type = str, 
                                      default = '0', 
                                      help = tools_id_helper)  
        
        # Tools combination command line argument:
        self.cmdLineArgs.add_argument('--comb', 
                                      type = str, 
                                      default = '[AND=]', 
                                      help = consensus_helper)  
        
        # Prediction analysis command line argument:
        self.cmdLineArgs.add_argument('--pred', 
                                      type = str, 
                                      default = 'no', 
                                      help = prediction_helper)
    
        # Enrichment analysis command line argument:
        self.cmdLineArgs.add_argument('--enrich', 
                                      type = str, 
                                      default = 'no', 
                                      help = enrichment_helper) 
        
        # Enrichment analysis type command line argument:
        self.cmdLineArgs.add_argument('--enrichtype', 
                                      type = str, 
                                      default = 'go', 
                                      help = enrichment_type_helper)        
        
        # Consensus command line argument:
        self.cmdLineArgs.add_argument('--cons', 
                                      type = int, 
                                      default = 1, 
                                      help = min_consensus_helper)  
        
        # Statistical test command line argument:
        self.cmdLineArgs.add_argument('--pvalue', 
                                      type = float, 
                                      default = 0.05, 
                                      help = pvalue_helper)  
        
        # Statistical test command line argument:
        self.cmdLineArgs.add_argument('--pvalue_type', 
                                      type = str, 
                                      default = 'raw',
                                      help = pvalue_type_helper)        
     
        # Set all command line argument into args
        args = self.cmdLineArgs.parse_args()      
        
        self.consensus = args.comb
        self.prediction = args.pred
        self.enrichment = args.enrich
        self.enrich_type = [ e for e in args.enrichtype.split(',') 
                             if e in ['go', 'kegg', 'reactome']]
        self.no_cpu_cores = args.cores
        self.min_consensus = args.cons
        self.pvalue = args.pvalue
        self.pvalue_type = args.pvalue_type
        self.tools = [ int(t.strip()) 
                       for t in args.tools.split(',') 
                       if isinstance(int(t), int) and int(t) in [0,1,2,3,4] ]  
        
        if self.prediction not in ['yes', 'no']:
            self.ePrintClose("Error." + prediction_helper)
        
        if self.enrichment not in ['yes', 'no']:
            self.ePrintClose("Error." + enrichment_helper)
                
        if not isinstance(self.no_cpu_cores, int): 
            self.ePrintClose("Error." + cpu_cores_helper)
         
        if not isinstance(self.min_consensus, int): 
            self.ePrintClose("Error." + min_consensus_helper)  
            
        if not self.tools:
            self.ePrintClose("Error." + tools_helper)
            
        if not self.enrich_type:
            self.ePrintClose("Error." + enrichment_type_helper)
                        
        if not isinstance(self.pvalue, float): 
            self.ePrintClose("Error." + pvalue_helper)
        
        if self.pvalue_type not in ['raw', 'adj']:
            self.ePrintClose("Error." + pvalue_type_helper)
            
        no_system_cores = self.available_cpu_count()
        
        # Checks the number of cores available
        if no_system_cores < self.no_cpu_cores:
            self.no_cpu_cores = no_system_cores
            
            cores_warning = '''
Warning. It has been specified a number of CPU cores greater than the available 
ones. The execution will use %d cores.\n\nNo. of cores used: %d
''' % (self.no_cpu_cores, self.no_cpu_cores )
                
                
    def ePrint(self, message):
        print >> sys.stderr, message
        
        
    def ePrintClose(self, message):
        print >> sys.stderr, message  
        sys.exit(1)    
        
        
    def available_cpu_count(self):
        """ Number of available virtual or physical CPUs on this system, i.e.
        user/real as output by time(1) when called with an optimally scaling
        userspace-only program"""
    
        # https://github.com/giampaolo/psutil
        try:
            import psutil
            return psutil.cpu_count()   # psutil.NUM_CPUS on old versions
        except (ImportError, AttributeError):
            pass
    
        # POSIX
        try:
            res = int(os.sysconf('SC_NPROCESSORS_ONLN'))
    
            if res > 0:
                return res
        except (AttributeError, ValueError):
            pass
    
        # Windows
        try:
            res = int(os.environ['NUMBER_OF_PROCESSORS'])
    
            if res > 0:
                return res
        except (KeyError, ValueError):
            pass
    
        return 1
            
    
