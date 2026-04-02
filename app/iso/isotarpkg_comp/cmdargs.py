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
        # Define the AgumentParser object
        self.cmdLineArgs = argparse.ArgumentParser(description = 
                                                       "Run run.py")
        
        cpu_cores_helper = '''Specify the right number of CPU cores to be used'''    
        list1_helper = '''Specify the first list of genes'''
        list2_helper = '''Specify the second list of genes'''
        pvalue_helper = '''Specify a pvalue'''
        pvalue_type_helper = '''Specify the pvalue type [raw, adj]'''
	enrichment_type_helper = '''Specify which type of enrichment 
analysis, as list, you would like to perform [go, kegg, reactome]'''	
        
        # No. of cores command line argument
        self.cmdLineArgs.add_argument('--cores', 
                                      type = int, 
                                      default = 1,
                                      help = cpu_cores_helper) 
        
        # First targets list
        self.cmdLineArgs.add_argument('--list1', 
                                      type = str, 
                                      required = "TRUE", 
                                      help = list1_helper)
        
        # Second targets list
        self.cmdLineArgs.add_argument('--list2', 
                                      type = str,
                                      default = '', 
                                      help = list2_helper) 
        
        # First targets list name
        self.cmdLineArgs.add_argument('--name1', 
                                      type = str, 
                                      required = "TRUE", 
                                      help = list1_helper)
        
        # Second targets list name
        self.cmdLineArgs.add_argument('--name2', 
                                      type = str,
                                      default = '', 
                                      help = list2_helper)        
        
        # P-value
        self.cmdLineArgs.add_argument('--pvalue', 
                                      type = float, 
                                      default = 0.05, 
                                      help = pvalue_helper)  
        
        # Statistical test
        self.cmdLineArgs.add_argument('--pvalue_type', 
                                      type = str, 
                                      default = 'raw', 
                                      help = pvalue_type_helper)  
    
	# Enrichment analysis type command line argument:
	self.cmdLineArgs.add_argument('--enrichtype', 
	                              type = str, 
	                              default = 'go', 
	                              help = enrichment_type_helper)
	
     
        # Set all command line argument into args
        args = self.cmdLineArgs.parse_args()      
        
        self.no_cpu_cores = args.cores
        self.pvalue = args.pvalue
        self.pvalue_type = args.pvalue_type
        self.targets_list1 = args.list1
        self.targets_list2 = args.list2
        self.name_list1 = args.name1
        self.name_list2 = args.name2    
	self.enrich_type = [ e for e in args.enrichtype.split(',') 
	                     if e in ['go', 'kegg', 'reactome']]	
                
        if not isinstance(self.no_cpu_cores, int): 
            self.ePrintClose("Error." + cpu_cores_helper)
            
        if self.targets_list1 == '':
            self.ePrintClose("Error." + list1_helper)
         
        if not isinstance(self.pvalue, float): 
            self.ePrintClose("Error." + pvalue_helper)
        
        if self.pvalue_type not in ['raw', 'adj']:
            self.ePrintClose("Error." + pvalue_type_helper)
	    
	if not self.enrich_type:
            self.ePrintClose("Error." + enrichment_type_helper)
            
        #no_system_cores = psutil.cpu_count()
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
		
    
