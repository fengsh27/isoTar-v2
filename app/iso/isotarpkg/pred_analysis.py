# -*- coding: utf-8 -*-

''' Libraries '''
import re
import argparse
import os, errno
import os.path
import sys
import psutil
import target as tr
import mirna as mr
import itutils as ut
import multiprocessing as mp
import json
import csv

class PredictionAnalysis():
        
    def __init__(self, mirna_path, target_path, output_dir, 
                 pre_mature_path, cmd_l, pre_mat_acc_path, nm_to_geneid_path):
        
        self.mirna_path = mirna_path
        self.target_path = target_path
        self.output_dir = output_dir
        self.pre_mature_path = pre_mature_path
	self.pre_mat_acc_path = pre_mat_acc_path
	self.nm_to_geneid_path = nm_to_geneid_path
        
        if not os.path.isfile(self.mirna_path):
            self.ePrintClose("Error. No mirnas.fa file found.")
        
        if not os.path.isfile(self.target_path):
            self.ePrintClose("Error. No 3utr.fa file found.")
        
        if not os.path.isdir(self.output_dir):
            self.ePrintClose("Error. No output directory found.")
            
        if not os.path.isfile(self.pre_mature_path):
            self.ePrintClose("Error. No pre_mature_path.json file found.")
            
        if not os.path.isfile(self.pre_mat_acc_path):
	    self.ePrintClose("Error. No pre_mat_acc.json file found.")
	    
	if not os.path.isfile(self.nm_to_geneid_path):
	    self.ePrintClose("Error. No nm_to_geneid.tab file found.")
        
        # //////////////////////////////////////////////////////////////////////
        # //////////////////////////////////////////////////////////////////////
        # miRNAs processing
        # //////////////////////////////////////////////////////////////////////
        # //////////////////////////////////////////////////////////////////////    
        
        # miRNA object
        self.mirna_obj = mr.MiRNAsProcessing(self.pre_mature_path,
	                                     self.pre_mat_acc_path)
        
        # miRNA loading
        self.mirna_obj.loadMiRNAs(self.mirna_path) 
        
        # //////////////////////////////////////////////////////////////////////
        # //////////////////////////////////////////////////////////////////////
        # Targets processing
        # //////////////////////////////////////////////////////////////////////
        # //////////////////////////////////////////////////////////////////////
        
        # Targets
        self.target_obj = tr.TargetsProcessing(target_path, 
                                               cmd_l.no_cpu_cores)
        
        # Targets processing
        self.target_obj.getTargets()        
        
        # //////////////////////////////////////////////////////////////////////
        # //////////////////////////////////////////////////////////////////////
        # Tools checking
        # //////////////////////////////////////////////////////////////////////
        # //////////////////////////////////////////////////////////////////////    
        
        # isoTar available tools list
        self.available_tools = [ 'miRanda', 
                                 'miRmap', 
                                 'TargetScan',
                                 'RNAhybrid',
                                 'PITA']
        
        self.selected_tools = []
        
        for i in xrange(0, len(self.available_tools)):
            if i in cmd_l.tools:
                tool = self.available_tools[i]
                if tool not in self.selected_tools:
                    self.selected_tools.append(tool)
		    
        
        # Split the 3utr file into N files (based on the no. of cores used)
        if any(item in ['RNAhybrid', 'PITA', 'miRmap', 'miRanda'] 
               for item in self.selected_tools):        
            self.target_obj.prepareTargets()    
            
        if 'TargetScan' in self.selected_tools:  
            
            ts_fam = '/opt/TargetScan/Datasets/miR_Family_Info.txt'
            ts_est_to_refseq = '/app/iso/static/public/resources/enst_to_refseq.txt'
                        
	    ts_3utr = '/opt/TargetScan/Datasets/3utr'
	    ts_bls_bins = '/opt/TargetScan/Datasets/bln_bins'  
            
            if not os.path.isfile(ts_fam):
                self.ePrintClose('Error. No TargetScan files found')
                
            if not os.path.isfile(ts_est_to_refseq):
                self.ePrintClose('Error. No Enst map file found')

            if not os.path.isdir(ts_3utr):
                self.ePrintClose('Error. No TargetScan 3utr directory found')
                
            if not os.path.isdir(ts_bls_bins):
                self.ePrintClose('Error. No TargetScan bls bins directory found')
                
            self.mirna_obj.getMirFamilyInfo(ts_fam)
            self.target_obj.targetsEnstToRefseqMap(ts_est_to_refseq)
            self.target_obj.prepareTargetscanTargets(ts_3utr, ts_bls_bins)

        self.intersection_tools = []
        intersection_tools = []
    
        # Tools intersection combination parsing
        matchObj = re.match(r'.*\[AND=([0-4,]+)\].*', 
                            cmd_l.consensus, 
                            re.M|re.I)
        if matchObj:
            intersection_tools = [ int(t.strip()) 
	                           for t in matchObj.group(1).split(',') 
	                           if isinstance(int(t), int) ]
	    
	self.intersection_tools = [ self.available_tools[el] for el in intersection_tools ]         
        
        # //////////////////////////////////////////////////////////////////////
        # Output directory setting
        # ////////////////////////////////////////////////////////////////////// 
        
        # The object provides a set of facilities to create and remove the output
        # directory, as well as functions to run bash commands and multi processes
        utl = ut.ITUtils(cmd_l.no_cpu_cores)
        
        # Delete an existing output direcotry
        utl.removeDirectoryContent('{0}/prediction'.format(self.output_dir))
        utl.removeDirectoryContent('{0}/enrichment'.format(self.output_dir))
	
        # Create a new empty output directory
        #utl.createDirectory(self.output_dir)
    
        # Get the targeting tools list
        self.tool_details = utl.getTargetingToolsList()
        
        # Workers queue
        tmp_queue = mp.Queue()
        pred_queue = mp.Queue()
        
        # miRNAs output directories creation
        self.mirna_obj.setMiRNAsOutputDirectories(self.output_dir)
                
        # Performe the targeting analysis for each selected tool
        results_list, final_results_list = utl.executeTargetingAnalysis(
            self.tool_details, self.selected_tools, 1, self.mirna_obj, 
            self.target_obj, tmp_queue, pred_queue)   
           
        # Remove all temporary splitted 3utr targets files
        self.target_obj.deleteTargetsSupportFile()           
        
        # Create a supporting file for the isoTar Web UI
	self.prepare_prediction_synthesis(cmd_l)
	
	
    def ePrint(self, message):
        print >> sys.stderr, message
        
    def ePrintClose(self, message):
        print >> sys.stderr, message  
        sys.exit(1)   
    
    def get_nm_to_geneid_map(self):
	nm_to_geneid_map = {}
	if os.path.isfile(self.nm_to_geneid_path):
	    with open(self.nm_to_geneid_path, 'r') as f:
		handler = csv.reader(f, delimiter='\t')
		next(handler, None)
		for line in handler:	
		    if len(line) == 3:
			nm_to_geneid_map[line[0]] = {
			    'accession_number': line[0],
			    'gene_id': line[1],
			    'gene_name': line[2]
			}
			
	return nm_to_geneid_map	    
    
    def prepare_prediction_synthesis(self, cmd_l):
        
        nm_to_geneid_map = self.get_nm_to_geneid_map()
	
        output_path = '%s/prediction' % self.output_dir
	syn = {}
	
	if os.path.isdir(output_path):
	    # Get all miRNA(s) prediction directories
	    predictions = next(os.walk(output_path))[1]
	    for prediction in predictions:
		mirna_json_path = '%s/%s/mirna_check.json' % (
		    output_path, 
		    prediction)  
		if os.path.isfile(mirna_json_path):
		    with open(mirna_json_path, 'r') as f:
			mirna_info = json.load(f)			
			if mirna_info and \
			   'identifier' in mirna_info and \
			   'output_path' in mirna_info:   
			    	    
			    mirna_id = mirna_info['identifier']
			    mirna_acc = mirna_info['acc']
			    
			    pre_id = mirna_info['pre_mirna_id']
			    pre_acc = mirna_info['acc']		
			    
			    if mirna_id not in syn:
				
				p, p_map = self.get_output_data_per_tool(
				    mirna_info['output_path'], 
				    cmd_l.min_consensus,
				    nm_to_geneid_map)
				
				syn[mirna_id] = {
				    'mirna_id': mirna_id,
				    'mirna_acc': mirna_acc,
				    'pre_acc': pre_acc,
				    'pre_mirna_id': mirna_info['pre_mirna_id'],
				    'variations': {
				        'substitution': [],
				        'isomir': []
				    },
				    'predictions': p,
				    'support': p_map
				}
				
				# Nucleotide substitution
				for s in mirna_info['variations']['substitution']:
				    if os.path.isdir(s['output_path']):
					
					p, p_map = self.get_output_data_per_tool(
					    s['output_path'], 
					    cmd_l.min_consensus,
					    nm_to_geneid_map)
					
					d = {
					    'position': s['position'],
					    'ref': s['ref'],
					    'mod': s['mod'],
					    'predictions': p,
					    'support': p_map					    
					} 
					syn[mirna_id]['variations']['substitution'].append(d)
				
				# isomiR
				for i in mirna_info['variations']['isomir']:
				    if os.path.isdir(i['output_path']):
					
					p, p_map = self.get_output_data_per_tool(
					    i['output_path'], 
					    cmd_l.min_consensus,
					    nm_to_geneid_map)
					
					d = {
					    '3p': i['3p'],
					    '5p': i['5p'],
					    'predictions': p,
					    'support': p_map
					} 
					syn[mirna_id]['variations']['isomir'].append(d)
				
	pred_json_path = '%s/prediction.json' % output_path	    
	with open(pred_json_path, 'w') as f:
	    f.write(json.dumps(syn, indent=4))	

    def get_output_data_per_tool(self, output_path, consensun, nm_to_geneid_map):
	
	output = {
	    'targets': [],
	    'miRanda': [],
	    'miRmap': [],
	    'TargetScan': [],
	    'RNAhybrid': [],
	    'PITA': []
	}
	
	output_map = {
	    'targets': [],
	    'miRanda': [],
	    'miRmap': [],
	    'TargetScan': [],
	    'RNAhybrid': [],
	    'PITA': []
	}	
	
	tools = ['miRanda', 'miRmap', 'TargetScan', 'RNAhybrid', 'PITA']
	
	path = '%s/output.tsv' % output_path
	
	if os.path.isfile(path):
	    with open(path, 'r') as f:
		handler = csv.reader(f, delimiter='\t')
		headers = handler.next()
		columns = {}
		for header in headers:
		    columns[header] = []
		
		for line in handler:
		    for h, v in zip(headers, line):
			columns[h].append(v)

		# Get the last occurrence IDX of the specified consensun
		if 'Sum' in columns:
		    # Get the latest row according to the specified consensus
		    last_target_idx = ''.join(columns['Sum']).rfind(str(consensun))
		    if last_target_idx != -1:
			# Create the global targets list
			output['targets'] = columns['Target'][0:last_target_idx+1]
			output_map['targets'] = [ nm_to_geneid_map[el] for el in output['targets'] if el in nm_to_geneid_map ]
			
			tools_combination_targets_idx = []
			
			# Check which tool must appear in our targets
			for tool in self.intersection_tools:
			    # Get all indexex in which the tool predicted a specific target
			    l1 = [m.start() for m in re.finditer('1', ''.join(columns[tool][0:last_target_idx+1]))]
			    # Map each index to the corresponding target
			    l2 = map(columns['Target'].__getitem__, l1)
			    tools_combination_targets_idx.append(l2)
			    
			# Update the global targets list
			if len(self.intersection_tools) > 0:
			    # Merge all lists into one
			    output['targets'] = list(set.intersection(*map(set, tools_combination_targets_idx)))
			    output_map['targets'] = [ nm_to_geneid_map[el] for el in output['targets'] if el in nm_to_geneid_map ]
			
			for tool in tools:
			    if tool in columns:
				# Get all indexex in which the tool predicted a specific target
				l = [m.start() for m in re.finditer('1', ''.join(columns[tool][0:last_target_idx+1]))]
				# Map each index to the corresponding target
				output[tool] = map(columns['Target'].__getitem__, l)
				output_map[tool] = [ nm_to_geneid_map[el] for el in output[tool] if el in nm_to_geneid_map ]
				
	return output, output_map