# -*- coding: utf-8 -*-

''' Libraries '''
import mimerender
import re
import datetime
from bson.json_util import dumps
import argparse
import os, errno
import os.path
import json
import sys
from glob import glob
from itertools import chain
import subprocess
import time
import shutil
import itutils as ut
import multiprocessing as mp
import math
from operator import itemgetter
import csv
from goatools.base import download_go_basic_obo
from goatools.base import download_ncbi_associations
from goatools.obo_parser import GODag
from goatools.associations import read_ncbi_gene2go
from goatools.test_data.genes_NCBI_9606_ProteinCoding import GeneID2nt as GeneID2nt_hsa
from goatools.go_enrichment import GOEnrichmentStudy

class GO():
        
    def __init__(self, nm_to_geneid_path, cmd_l, go_path, output_path):
        
        if not os.path.isfile(nm_to_geneid_path):
            self.ePrintClose("Error. No nm_to_geneid_path.tab file found.")
	
	if not os.path.isfile(go_path):
            self.ePrintClose("Error. No go.json file found.")
	
	if not os.path.isdir(output_path):
            self.ePrintClose("Error. No enrichment_comp directory found.")
	    
	self.output_enrichment_dir = output_path
        self.debug = False          
        self.nm_to_geneid_path = nm_to_geneid_path
        self.refseq_gene_name_map = {}
	self.gene_name_gene_id_map = {}
	self.pvalue = cmd_l.pvalue
	self.updated_dot_files_path_list = []
	self.go_path = go_path
	self.go_info = {}
	self.pvalue_type = cmd_l.pvalue_type
	self.input_targets_list1 = [ el.strip() for el in cmd_l.targets_list1.replace('"', '').replace('\'', '').split(',') if el.strip() != '' ]  
	self.input_targets_list2 = [ el.strip() for el in cmd_l.targets_list2.replace('"', '').replace('\'', '').split(',') if el.strip() != '' ]	
	self.discarded_genes_in_list1 = []
	self.discarded_genes_in_list2 = []
	self.targets_list1 = []  
	self.targets_list2 = []
	self.name_list1 = re.sub(r'\s+', '_', re.sub(r'[^0-9a-zA-Z_\s\-]+', '', cmd_l.name_list1))  
	self.name_list2 = re.sub(r'\s+', '_', re.sub(r'[^0-9a-zA-Z_\s\-]+', '', cmd_l.name_list2)) 	
	self.support = {}
	self.give_me_my_dot_file = 5
	self.geneid_to_genename = {}
		
	with open(self.go_path, 'r') as f:
	    self.go_info = json.load(f)
	
	
    def ePrint(self, message):
        print >> sys.stderr, message
        
        
    def ePrintClose(self, message):
        print >> sys.stderr, message  
        sys.exit(1)
        
        
    def setDebugging(self, debug):
        self.debug = debug    
    

    def executeCommandLine(self, command):
        retcode = 0
        try:
            retcode = subprocess.check_output(command, shell = True)
            
            if self.debug:
                if retcode < 0:
                    print >>sys.stderr, "Child was terminated by signal", -retcode
                else:
                    print >>sys.stderr, "Child returned", retcode
                    pass
	except subprocess.CalledProcessError as e:
	    if self.debug:
		print >>sys.stderr, "Execution failed:", e.returncode	 
	    pass         
        except OSError as e:
            if self.debug:
                print >>sys.stderr, "Execution failed:", e
	    pass
            
        return retcode  
    
                        
    def getGeneIdMap(self): 
        
        if os.path.isfile(self.nm_to_geneid_path):
            with open(self.nm_to_geneid_path, 'r') as f:
                handler = csv.reader(f, delimiter='\t')
                # Skip the header
                next(handler, None)
                for line in handler:
                    # refseq geneid genename
                    if len(line) == 3 and line[0] != '' and line[1] != '' and \
                       line[2] != '':
                        refseq = line[0]
                        gene_id = line[1]
                        gene_name = line[2]
                        
                        if refseq not in self.refseq_gene_name_map:
                            self.refseq_gene_name_map[refseq] = {
                                'gene_name': gene_name,
                                'gene_id': gene_id
                            } 
			
			if gene_name not in self.gene_name_gene_id_map:
                            self.gene_name_gene_id_map[gene_name] = gene_id
			
			if gene_id not in self.geneid_to_genename:
			    self.geneid_to_genename[gene_id] = gene_name
			    
			    
    def checkInputGenesList(self): 
		
	# Check the first genes list
        for gene in self.input_targets_list1:
	    gene_name = ''
	    if gene in self.refseq_gene_name_map:
		g_name = self.refseq_gene_name_map[gene]['gene_name']
		if g_name in self.gene_name_gene_id_map:
		    gene_name = g_name
	    elif gene in self.gene_name_gene_id_map:
		gene_name = gene
	    
	    # Check if the gene name has been set
	    if gene_name == '':
		# Discard it!
		self.discarded_genes_in_list1.append(gene)
	    elif gene_name not in self.targets_list1:
		self.targets_list1.append(gene_name)
			
	# Check the second genes list
        for gene in self.input_targets_list2:
	    gene_name = ''
	    if gene in self.refseq_gene_name_map:
		g_name = self.refseq_gene_name_map[gene]['gene_name']
		if g_name in self.gene_name_gene_id_map:
		    gene_name = g_name
	    elif gene in self.gene_name_gene_id_map:
		gene_name = gene
	    
	    # Check if the gene name has been set
	    if gene_name == '':
		# Discard it!
		self.discarded_genes_in_list2.append(gene)
	    elif gene_name not in self.targets_list2:
		self.targets_list2.append(gene_name)
	
	self.support = {
	    'name_list1': self.name_list1,
	    'name_list2': self.name_list2,	    
	    'targets_list1': self.targets_list1,
	    'targets_list2': self.targets_list2,
	    'targets_list1_discarded': self.discarded_genes_in_list1,
	    'targets_list2_discarded': self.discarded_genes_in_list2	    
	}
	
	# Restrictions
	if self.targets_list1 == self.targets_list2:
	    self.targets_list2 = []
	    
	if len(self.targets_list2) == 0:
	    self.name_list2 = ''
	
	if self.name_list2 == '':
	    self.targets_list2 = []	
	    
	
	                            
                                    
    def executeEnrichmentAnalysis(self, no_cpu_cores, queue, obo_dag, 
                                  geneid2gos_human, utils):
        # Create the gene id map
        self.getGeneIdMap()
	
	# Check genes input lists
	self.checkInputGenesList()
        
        # Create 
        enrich_dirs_per_core = [ [] for i in xrange(0, no_cpu_cores) ]
	enrich_dirs = []
	
	# Dir labels
	comparison_dirs = ['', '', '']
	
	'''
	////////////////////////////////////////////////////////////////////////
	First targets list
	////////////////////////////////////////////////////////////////////////
	'''
	if len(self.targets_list1) > 0:
	    dir_path = '{0}/{1}'.format(self.output_enrichment_dir, 
	                                self.name_list1)
	    enrich_dirs.append({
	        'enrichment_dir': self.name_list1,
	        'targets': self.targets_list1,
	        'enrichment_dir_path': dir_path,
	        'label': 'list1',
	        'type': 'single'
	    })
	    utils.createDirectory(dir_path)
	    self.support['enrichment_list1_dir_path'] = dir_path

	'''
	////////////////////////////////////////////////////////////////////////
	Second targets list
	////////////////////////////////////////////////////////////////////////
	'''    
	if len(self.targets_list2) > 0:
	    dir_path = '{0}/{1}'.format(self.output_enrichment_dir, 
	                                self.name_list2)	    
	    enrich_dirs.append({
	        'enrichment_dir': self.name_list2,
	        'targets': self.targets_list2,
	        'enrichment_dir_path': dir_path,
	        'label': 'list2',
	        'type': 'single'	        
	    })
	    utils.createDirectory(dir_path)
	    self.support['enrichment_list2_dir_path'] = dir_path
	    
	    '''
	    ////////////////////////////////////////////////////////////////////
	    Targets beloging to the first list only
	    ////////////////////////////////////////////////////////////////////
	    '''	    
	    targets_list1_only = list(set(self.targets_list1) - 
	                              set(self.targets_list2))
	    
	    # Targets list1 only dir path
	    dir_path = '{0}/comparison/{1}_only'.format(
	        self.output_enrichment_dir, self.name_list1)
	    utils.createDirectory(dir_path)
	    self.support['enrichment_comp_list1_only_dir_path'] = dir_path
	    self.support['enrichment_comp_dir_path'] = '{0}/comparison'.format(
	        self.output_enrichment_dir)
	    
	    if len(targets_list1_only) > 0:
		enrich_dirs.append({
		    'enrichment_dir': '{0}_only'.format(self.name_list1),
		    'targets': targets_list1_only,
		    'enrichment_dir_path': dir_path,
		    'label': 'list1_only',
		    'type': 'comparison'
		})
		comparison_dirs[0] = '{0}_only'.format(self.name_list1)
		
	    
	    '''
	    ////////////////////////////////////////////////////////////////////
	    Targets beloging to the second list only
	    ////////////////////////////////////////////////////////////////////
	    '''	
	    targets_list2_only = list(set(self.targets_list2) - 
	                              set(self.targets_list1))
	    
	    # Targets list2 only dir path
	    dir_path = '{0}/comparison/{1}_only'.format(
	        self.output_enrichment_dir, self.name_list2)	    
	    utils.createDirectory(dir_path)
	    self.support['enrichment_comp_list2_only_dir_path'] = dir_path
	    
	    if len(targets_list2_only) > 0:
		enrich_dirs.append({
		    'enrichment_dir': '{0}_only'.format(self.name_list2),
		    'targets': targets_list2_only,
		    'enrichment_dir_path': dir_path,
		    'label': 'list2_only',
		    'type': 'comparison'
		})
		comparison_dirs[1] = '{0}_only'.format(self.name_list2)
	    
	    '''
	    ////////////////////////////////////////////////////////////////////
	    Targets beloging to the intersection between first and second list
	    ////////////////////////////////////////////////////////////////////
	    '''
	    # Intersection between genes list1 and list2
	    intersec_targets = list(set(self.targets_list1).intersection(
	        set(self.targets_list2)))	
	    
	    # Intersection path
	    dir_path = '{0}/comparison/intersection'.format(self.output_enrichment_dir)	  
	    utils.createDirectory(dir_path)
	    self.support['enrichment_comp_intersection_dir_path'] = dir_path
	    
	    if len(intersec_targets) > 0:
		enrich_dirs.append({
		    'enrichment_dir': 'intersection',
		    'targets': intersec_targets,
		    'enrichment_dir_path': dir_path,
		    'label': 'intersection',
		    'type': 'comparison'
		})
		comparison_dirs[2] = 'intersection'
	        
        # Prepare all enrichment files to be assigned to each available core
        for i in xrange(0, len(enrich_dirs)):
            idx = int(math.fmod(i, no_cpu_cores))
            # Add the i-th enrichment to the idx-th core
            enrich_dirs_per_core[idx].append(enrich_dirs[i])
        
        jobs = []
        
        # For each core
        for i in xrange(0, len(enrich_dirs_per_core)):
            if enrich_dirs_per_core[i]:
                p = mp.Process(target = self.worker, 
                               args=(i,
                                     enrich_dirs_per_core[i], # A pool of prediction(s) dir
                                     queue,
                                     obo_dag, 
		                     geneid2gos_human))               
                jobs.append(p)
                p.start() 
        
        # Wait worker(s) termination
        self.waitWorker(jobs, queue)  
	
	'''
	////////////////////////////////////////////////////////////////////////
	Comparison between list 1 and list 2
	////////////////////////////////////////////////////////////////////////
	'''
	
	if len(self.targets_list1) > 0 and len(self.targets_list2) > 0:
	    # GO terms Tree
	    trees_list = self.run_r_go_term_comparison_tree(
		'{0}/comparison'.format(self.output_enrichment_dir), comparison_dirs)
	    
	    for tree in trees_list:
		# Make the GO terms Tree
		self.create_go_terms_tree(tree['dot_path'], 
		                          tree['pdf_path'],
		                          tree['work_dir_path'])	
        
                              
    def worker(self, idx, enrichment_dirs, queue, obo_dag, geneid2gos_human):
        
        # For each enrichment dir
        for enrichment in enrichment_dirs:
	    
	    trees_list = []
            if len(enrichment['targets']) > 0:
                # Enrichment 
                self.computeEnrichment(enrichment['targets'], 
                                       enrichment['enrichment_dir_path'],
                                       obo_dag, 
		                       geneid2gos_human)
		
		if enrichment['type'] == 'single':
		    # GO terms Tree
		    trees_list = self.run_r_go_term_dot_file(
		        enrichment['enrichment_dir_path'],
		        enrichment['enrichment_dir'])
            
            queue.put(({ 
                'worker': idx, 
	        'trees_list': trees_list
            }, idx))
        
    
    def waitWorker(self, jobs, queue):
        
        trees_list = []
        # Collects all results into in a list                    
        for j in xrange(0, len(jobs)):
            # Get worker result
            res = queue.get()
            if res:
                worker_id = res[0]['worker']
		if len(res[0]['trees_list']) > 0:
		    trees_list += res[0]['trees_list']
                
        # Wait for all worker processes to finish
        for p in jobs:
            p.join() 
	
	for tree in trees_list:
	    # Make the GO terms Tree
	    self.create_go_terms_tree(tree['dot_path'], 
                                      tree['pdf_path'],
                                      tree['work_dir_path']) 
	    
	    
    def computeEnrichment(self, targets_list, path, obo_dag, geneid2gos_human):
        
	goea_obj = GOEnrichmentStudy(
            GeneID2nt_hsa.keys(), # List of human protein-coding genes
            geneid2gos_human, # geneid/GO associations
            obo_dag, # Ontologies
            propagate_counts = False,
            alpha = 0.05, # default significance cut-off
            methods = ['fdr_bh']) # defult multipletest correction method
	
        header = '\t'.join(['Go term', 'GO name', 'GO namespace', 'P-value', 
                            '-Log(P-value)', 'Adjusted P-value', 
                            '-Log(Adjusted P-value)', 'Ratio in study', 
                            'Ratio in population', 'Gene names'])
        header += '\n'    
	
	# ======================================================================
	# TSV
	# ======================================================================	
        # Enrichment cellular component file path
        enrich_cc_f_path = '%s/enrich_cc.tsv' % path
    
        # Enrichment biological process file path
        enrich_bp_f_path = '%s/enrich_bp.tsv' % path
       
        # Enrichment molecular function file path
        enrich_mf_f_path = '%s/enrich_mf.tsv' % path
	
	# ======================================================================
	# JSON
	# ======================================================================
	# Enrichment cellular component file path
	enrich_cc_f_json_path = '%s/enrich_cc.json' % path
    
	# Enrichment biological process file path
	enrich_bp_f_json_path = '%s/enrich_bp.json' % path
       
	# Enrichment molecular function file path
	enrich_mf_f_json_path = '%s/enrich_mf.json' % path
	
	# ======================================================================
	# JSON GO list
	# ======================================================================
	# Enrichment cellular component go list file path
	enrich_cc_f_go_json_path = '%s/enrich_cc_go.json' % path
    
	# Enrichment biological process go list file path
	enrich_bp_f_go_json_path = '%s/enrich_bp_go.json' % path
       
	# Enrichment molecular function go list file path
	enrich_mf_f_go_json_path = '%s/enrich_mf_go.json' % path	
	
	enrich_cc = []
	enrich_bp = []
	enrich_mf = []
	
	enrich_cc_json = {}
	enrich_bp_json = {}
	enrich_mf_json = {}
	
	enrich_cc_go_json = {}
	enrich_bp_go_json = {}
	enrich_mf_go_json = {}
                
        if targets_list:
            input_genes = {}    
            
            for gene_name in targets_list:
		gene_id = int(self.gene_name_gene_id_map[gene_name])
		input_genes[gene_id] = gene_name
            
            geneids_study = input_genes.keys()
            goea_results_all = goea_obj.run_study(geneids_study)  
	    	                
            for r in goea_results_all:
		gene_names_list = [ self.geneid_to_genename[str(gene_id)] for gene_id in list(r.study_items) if str(gene_id) in self.geneid_to_genename ] if len(list(r.study_items)) > 0 else []
		
		if gene_names_list:
		    if r.NS == 'BP':
			l = [
			    r.GO, 
			    r.name, 
			    'biological_process', 
			    '{0:5.2e}'.format(r.p_uncorrected), 
			    '{:f}'.format(math.fabs(math.log10(r.p_uncorrected))), 
			    '{0:5.2e}'.format(r.p_fdr_bh), 
			    '{:f}'.format(math.fabs(math.log10(r.p_fdr_bh))),
			    '{0}/{1}'.format(r.study_count, r.study_n), 
			    '{0}/{1}'.format(r.pop_count, r.pop_n),
			    '{0}'.format(';'.join(gene_names_list))
			]		    
			enrich_bp.append('\t'.join(l) + '\n') 
			enrich_bp_go_json[r.GO] = math.fabs(math.log10(r.p_uncorrected))
			
		    elif r.NS == 'MF':
			l = [
			    r.GO, 
			    r.name, 
			    'molecular_function', 
			    '{0:5.2e}'.format(r.p_uncorrected), 
			    '{:f}'.format(math.fabs(math.log10(r.p_uncorrected))), 
			    '{0:5.2e}'.format(r.p_fdr_bh), 
			    '{:f}'.format(math.fabs(math.log10(r.p_fdr_bh))),
			    '{0}/{1}'.format(r.study_count, r.study_n), 
			    '{0}/{1}'.format(r.pop_count, r.pop_n),
			    '{0}'.format(';'.join(gene_names_list))
			]		    
			enrich_mf.append('\t'.join(l) + '\n') 
			enrich_mf_go_json[r.GO] = math.fabs(math.log10(r.p_uncorrected))
			
		    elif r.NS == 'CC':
			l = [
			    r.GO, 
			    r.name, 
			    'cellular_component', 
			    '{0:5.2e}'.format(r.p_uncorrected), 
			    '{:f}'.format(math.fabs(math.log10(r.p_uncorrected))), 
			    '{0:5.2e}'.format(r.p_fdr_bh), 
			    '{:f}'.format(math.fabs(math.log10(r.p_fdr_bh))),
			    '{0}/{1}'.format(r.study_count, r.study_n), 
			    '{0}/{1}'.format(r.pop_count, r.pop_n),
			    '{0}'.format(';'.join(gene_names_list))
			]		    
			enrich_cc.append('\t'.join(l) + '\n')  
			enrich_cc_go_json[r.GO] = math.fabs(math.log10(r.p_uncorrected))
		    
		    # Use a specific pvalue type: raw or adj
		    pvalue = float(r.p_uncorrected)
		    if self.pvalue_type == 'adj':
			pvalue = float(r.p_fdr_bh)
			
		    # P-value filtering according to the user
		    if pvalue < self.pvalue:
			if r.NS == 'BP':
			    d = {
				'go_name': r.name, 
				'go_namespace': 'biological_process', 
				'pvalue': '{0:5.2e}'.format(r.p_uncorrected), 
				'-Log(pvalue)': '{:f}'.format(math.fabs(math.log10(r.p_uncorrected))), 
				'adj_pvalue': '{0:5.2e}'.format(r.p_fdr_bh), 
				'-Log(adj_pvalue)': '{:f}'.format(math.fabs(math.log10(r.p_fdr_bh))),
				'ratio_in_study': '{0}/{1}'.format(r.study_count, r.study_n), 
				'ratio_in_population': '{0}/{1}'.format(r.pop_count, r.pop_n),
				'targets_list': gene_names_list
			    }
			    enrich_bp_json[r.GO] = d
			
			elif r.NS == 'MF':
			    d = {
				'go_name': r.name, 
				'go_namespace': 'molecular_function', 
				'pvalue': '{0:5.2e}'.format(r.p_uncorrected), 
				'-Log(pvalue)': '{:f}'.format(math.fabs(math.log10(r.p_uncorrected))), 
				'adj_pvalue': '{0:5.2e}'.format(r.p_fdr_bh), 
				'-Log(adj_pvalue)': '{:f}'.format(math.fabs(math.log10(r.p_fdr_bh))),
				'ratio_in_study': '{0}/{1}'.format(r.study_count, r.study_n), 
				'ratio_in_population': '{0}/{1}'.format(r.pop_count, r.pop_n),
				'targets_list': gene_names_list
			    }			
			    enrich_mf_json[r.GO] = d
			    
			elif r.NS == 'CC':
			    d = {
				'go_name': r.name, 
				'go_namespace': 'cellular_function', 
				'pvalue': '{0:5.2e}'.format(r.p_uncorrected), 
				'-Log(pvalue)': '{:f}'.format(math.fabs(math.log10(r.p_uncorrected))), 
				'adj_pvalue': '{0:5.2e}'.format(r.p_fdr_bh), 
				'-Log(adj_pvalue)': '{:f}'.format(math.fabs(math.log10(r.p_fdr_bh))),
				'ratio_in_study': '{0}/{1}'.format(r.study_count, r.study_n), 
				'ratio_in_population': '{0}/{1}'.format(r.pop_count, r.pop_n),
				'targets_list': gene_names_list
			    }			
			    enrich_cc_json[r.GO] = d
            
            # ==================================================================
	    # TSV
	    # ==================================================================
            with open(enrich_cc_f_path, 'w') as f:
                f.write(header)
                f.write(''.join(enrich_cc))
              
            with open(enrich_bp_f_path, 'w') as f:
                f.write(header)
                f.write(''.join(enrich_bp))
              
            with open(enrich_mf_f_path, 'w') as f:
                f.write(header)
                f.write(''.join(enrich_mf))
	    
	    # ==================================================================
	    # JSON
	    # ==================================================================
	    with open(enrich_cc_f_json_path, 'w') as f:
		f.write(json.dumps(enrich_cc_json, indent=4)) 
	      
	    with open(enrich_bp_f_json_path, 'w') as f:
		f.write(json.dumps(enrich_bp_json, indent=4)) 
	      
	    with open(enrich_mf_f_json_path, 'w') as f:
		f.write(json.dumps(enrich_mf_json, indent=4)) 
		
	    # ==================================================================
	    # JSON GO list
	    # ==================================================================
	    with open(enrich_cc_f_go_json_path, 'w') as f:
		f.write(json.dumps(enrich_cc_go_json, indent=4)) 
	      
	    with open(enrich_bp_f_go_json_path, 'w') as f:
		f.write(json.dumps(enrich_bp_go_json, indent=4)) 
	      
	    with open(enrich_mf_f_go_json_path, 'w') as f:
		f.write(json.dumps(enrich_mf_go_json, indent=4))	    
        else:
	    # ==================================================================
	    # TSV
	    # ==================================================================	    
            with open(enrich_cc_f_path, 'w') as f:
                f.write(header)
              
            with open(enrich_bp_f_path, 'w') as f:
                f.write(header)
              
            with open(enrich_mf_f_path, 'w') as f:
                f.write(header)
	    
	    # ==================================================================
	    # JSON
	    # ==================================================================
	    with open(enrich_cc_f_json_path, 'w') as f:
		f.write(json.dumps(enrich_cc_json, indent=4)) 
	      
	    with open(enrich_bp_f_json_path, 'w') as f:
		f.write(json.dumps(enrich_bp_json, indent=4)) 
	      
	    with open(enrich_mf_f_json_path, 'w') as f:
		f.write(json.dumps(enrich_mf_json, indent=4))
		
	    # ==================================================================
	    # JSON GO list
	    # ==================================================================
	    with open(enrich_cc_f_go_json_path, 'w') as f:
		f.write(json.dumps(enrich_cc_go_json, indent=4)) 
	      
	    with open(enrich_bp_f_go_json_path, 'w') as f:
		f.write(json.dumps(enrich_bp_go_json, indent=4)) 
	      
	    with open(enrich_mf_f_go_json_path, 'w') as f:
		f.write(json.dumps(enrich_mf_go_json, indent=4))
    
    
    def create_go_terms_dot_file(self, path, go_terms_list, work_dir_path):
	if len(go_terms_list) > 0:
	    command = '/usr/bin/Rscript /app/iso/isotarpkg_comp/create_go_terms_tree.R --input {0} --go {1} --work {2}'.format(
	        path, ','.join(go_terms_list), work_dir_path)
	    i = 0
	    while i < self.give_me_my_dot_file:
		self.executeCommandLine(command) 
		
		if os.path.isfile(path):
		    break
		time.sleep(1)
		i += 1	    


    def create_go_terms_tree(self, dot_path, pdf_path, work_dir_path):
	if os.path.isfile(dot_path):
	    command = '/usr/bin/Rscript /app/iso/isotarpkg_comp/create_updated_go_terms_tree.R --input {0} --output {1} --work {2}'.format(
	        dot_path, pdf_path, work_dir_path)
	    self.executeCommandLine(command)  
	    
	    
    def run_r_go_term_dot_file(self, path, label):
	
	go_class = ['bp', 'cc', 'mf']
	trees_list = []
	    
	# Extract from each JSON the GO terms list
        for go_c in go_class:
	    enrich_f_json_path = '{0}/enrich_{1}.json'.format(path, go_c)
	    if os.path.isfile(enrich_f_json_path):
		with open(enrich_f_json_path, 'r') as f:
		    d = json.load(f)
		    f_path = '{0}/go_tree_{1}.dot'.format(path, go_c)
		    
		    # Make the GO terms Tree
		    self.create_go_terms_dot_file(f_path, d.keys(), path)
		    
		    # GO terms color list
		    go_color_pvalue_list = { g: {
		        'color': '#1dc942',
		        'pvalue': d[g]['pvalue']
		    } for g in d.keys() }	
		    
		    input_dot_file = '{0}/go_tree_{1}.dot'.format(path, go_c)
		    output_dot_file = '{0}/go_tree_update_{1}.dot'.format(path, go_c)	
		    output_pdf_file = '{0}/{1}_go_tree_{2}.pdf'.format(
		        path, label, go_c)
		    	    
		    # Update the Dot file and create te GO term Tree
		    trees_list += self.update_go_tree_dot_file(input_dot_file, 
		                                               output_dot_file,
		                                               output_pdf_file,
		                                               go_color_pvalue_list,
		                                               path)
		    
	return trees_list
		    
		    
    def run_r_go_term_comparison_tree(self, comp_path, comparison_dirs):
	
	# Comparison directory path
	go_class = ['bp', 'cc', 'mf']
	trees_list = []
	
	colors = {
	    '001': '#ff6347', # List2: Red
	    '010': '#3cb371', # Intersection: Green
	    '011': '#9e8b5c', # Intersection & Variation:Brown
	    '100': '#1e90ff', # List1: Blue
	    '101': '#8f7aa3', # List1 & List2: Violet
	    '110': '#2da2b8', # List1 & Intersection: Dodger Blue
	    '111': '#b4b4b4' # All: Gray
	}
	    
	# Extract from each JSON the GO terms list
        for go_c in go_class:
	    go_list = []
	    
	    comp = {
	        'list1_only': [],
	        'list2_only': [],
	        'intersection': []
	    }
	    
	    comp_pvalue = {
	        'list1_only': {},
	        'list2_only': {},
	        'intersection': {}	        
	    }
	    
	    comp_dir = {
	        'list1_only': comparison_dirs[0],
	        'list2_only': comparison_dirs[1],
	        'intersection': comparison_dirs[2]
	    }
	    
	    # Retrieve GO terms from Wild Type only, Intersection, and Variation only
	    for c, l in comp.iteritems():
		if comp_dir[c] != '':
		    enrich_f_json_path = '{0}/{1}/enrich_{2}.json'.format(
			comp_path, comp_dir[c], go_c)
		    if os.path.isfile(enrich_f_json_path):
			with open(enrich_f_json_path, 'r') as f:
			    d = json.load(f)
			    comp[c] = d.keys()			
			    comp_pvalue[c] = { g: d[g]['pvalue'] for g in d.keys() }
	    
	    go_list = list(set(comp['list1_only']+comp['intersection']+comp['list2_only']))
	    
	    # Dot file path
	    f_path = '{0}/go_tree_{1}.dot'.format(comp_path, go_c)
	    
	    # Make the GO terms Tree
	    self.create_go_terms_dot_file(f_path, go_list, comp_path)
	    
	    go_color_pvalue_list = {}
	    
	    # For each GO term checks in which list the GO  term is
	    for g in go_list:
		b1 = '1' if g in comp['list1_only'] else '0'
		b2 = '1' if g in comp['intersection'] else '0'
		b3 = '1' if g in comp['list2_only'] else '0'
		
		label = '{0}{1}{2}'.format(b1, b2, b3)
		   
		if label in colors:
		    # GO terms color list
		    go_color_pvalue_list[g] = {
		        'color': colors[label],
		        'pvalue_list1_o': comp_pvalue['list1_only'][g] if g in comp_pvalue['list1_only'] else '',
		        'pvalue_inter': comp_pvalue['intersection'][g] if g in comp_pvalue['intersection'] else '',
		        'pvalue_list2_o': comp_pvalue['list2_only'][g] if g in comp_pvalue['list2_only'] else ''
		    }
		
	    input_dot_file = '{0}/go_tree_{1}.dot'.format(comp_path, go_c)
	    output_dot_file = '{0}/go_tree_update_{1}.dot'.format(comp_path, go_c)	
	    output_pdf_file = '{0}/go_tree_comp_{1}.pdf'.format(comp_path, go_c)
	    		    
	    # Update the Dot file and create te GO term Tree
	    trees_list += self.update_go_tree_dot_file(input_dot_file, 
	                                               output_dot_file,
	                                               output_pdf_file,
	                                               go_color_pvalue_list,
	                                               comp_path,
	                                               comparison_dirs)
	    
	return trees_list
	    
    
    def get_filtered_digraph(self, lines):
	
	u_lines = []
	if not lines:
	    return ['digraph test {',
	            '        graph [ratio=fill, size="12,8"];',
	            '        node [color=black,fillcolor=white,fontcolor=black,fontsize=10,label="\N",shape=plaintext,style=filled];',
	            '       edge [fontsize=10];',
	            '}']
	
	for i in range(len(lines)):
	    obj1 = re.match(r'^digraph.*$', lines[i], re.M| re.I)
	    obj2 = re.match(r'^\s+graph.*$', lines[i], re.M| re.I)
	    obj3 = re.match(r'^\s+node\s+.*$', lines[i], re.M| re.I)
	    obj4 = re.match(r'^\s+edge\s+.*$', lines[i], re.M| re.I)
	    obj5 = re.match(r'^\s+node[0-9]+\s+\[.*$', lines[i], re.M| re.I)
	    obj6 = re.match(r'^\s+node[0-9]+\s+\-\>\s+node[0-9]+\s+\[.*$', lines[i], re.M| re.I)
	    
	    # Check digraph
	    if obj1: 
		u_lines.append('digraph "GO terms Tree" {')
	    # Check graph
	    elif obj2:
		obj2a = re.match(r'^\s+.*\];$', lines[i], re.M| re.I)
		if not obj2a:
		    while not obj2a:
			i += 1
			obj2a = re.match(r'^\s+.*\];$', lines[i], re.M| re.I)
		u_lines.append('\tgraph [size="90,60", ratio=fill];')
	    # Check node
	    elif obj3:
		obj3a = re.match(r'^\s+.*\];$', lines[i], re.M| re.I)
		if not obj3a:
		    while not obj3a:
			i += 1
			obj3a = re.match(r'^\s+.*\];$', lines[i], re.M| re.I)
		u_lines.append('\tnode [label="\N", color=black, fillcolor=white, fontcolor=black, fontsize=8, shape=ellipse, style=filled];')
	    # Check edge
	    elif obj4:
		obj4a = re.match(r'^\s+.*\];$', lines[i], re.M| re.I)
		if not obj4a:
		    while not obj4a:
			i += 1
			obj4a = re.match(r'^\s+.*\];$', lines[i], re.M| re.I)		    
		u_lines.append('\tedge [fontsize=10];')
	    # Check nodeX or nodeX -> nodeY, with X = 1, 2,.., N
	    elif obj5 or obj6:
		line = ''
		obj5a = re.match(r'^\s+.*\];$', lines[i], re.M| re.I)
		if obj5a:
		    line = lines[i].strip()
		else:
		    while not obj5a:
			line += lines[i].strip()
			i += 1
			obj5a = re.match(r'^\s+.*\];$', lines[i], re.M| re.I)
		    line += lines[i].strip()
		u_lines.append('\t{}\n'.format(line))
	    else:
		u_lines.append(lines[i])
	
	return u_lines
    
    
    def update_go_tree_dot_file(self, input_dot_file, output_dot_file,
                                output_pdf_file, go_color_pvalue_list,
                                work_dir_path, labels = []):
	trees_list = []
	
	if os.path.isfile(input_dot_file) and len(go_color_pvalue_list) > 0:
	    with open(output_dot_file, 'w') as o_f:
		with open(input_dot_file, 'r') as i_f:
		    data = i_f.readlines()
		    
		    # Prepare all lines for the filtering
		    lines = self.get_filtered_digraph(data)
		    
		    # Dot file header
		    o_f.write('digraph "GO terms Tree" {\n')
		    o_f.write('\tgraph [size="90,60", ratio=fill];\n')
		    o_f.write('\tnode [label="\N", color=black, fillcolor=white, fontcolor=black, fontsize=8, shape=ellipse, style=filled];\n')
		    o_f.write('\tedge [fontsize=10];\n')
		    
		    # For each line into dot file
		    for i in xrange(4, len(lines)):
			line = ''
			# Check the node line
			matchObj1 = re.match(r'^\s+node[0-9]+\s+\[.+TD\>([a-z]+:[a-z0-9]+)\<br\/\>.+\];$', lines[i], re.M| re.I)
			if matchObj1: 
			    node_go_term = matchObj1.group(1)
			    node_idx = ''
			    matchObj2 = re.match(r'^\s+(node[0-9]+)\s+\[.+TD\>[a-z]+:[a-z0-9]+.+\];$', lines[i], re.M| re.I)
			    if matchObj2:
				node_idx = matchObj2.group(1)
			    '''
			    RamiGO bug: it does not escape special characters like &, ", <, > 
			    available within the GO term name, resulting in a break of the dot file.
			    Let's handle this!
			    
			    URL: https://www.graphviz.org/doc/info/shapes.html
			    
			    E.g.: GO term GO:0046920 has special character > in its name
			    
			    node9 [label=<<TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0"><TR><TD>GO:0046920<br/>alpha-(1->, color="#000000", fillcolor="#ffffff", fontcolor="#000000", 3=true];
			    fucosyl;
			    <br />;
			    ransferase;
			    activity;
			    </TD>;
			    </TR>;
			    </TABLE>;
			    node10
			    '''
			    if lines[i].find('</TD></TR></TABLE>>') == -1:
				'''
				\o/
				 |   Houston, we have a problem
				/ \
				'''
				# Escape all bad lines
				if i < (len(lines)-1):
				    i += 1
				    while not re.match(r'^\s+node[0-9]+\s+\[.+TD\>([a-z]+:[a-z0-9]+)\<br\/\>.+\];$', lines[i], re.M| re.I):
					i += 1
					
				    # For the next node
				    i -= 1
					
				if node_go_term in self.go_info:
				    go_name = self.go_info[node_go_term]['name']
				    # Replace bad characters
				    go_name = go_name.replace('&', '&amp;').replace('"', '&#34;').replace('<', '&lt;').replace('>', '&gt;')
				    
				    line = '\t{0} [label=<<TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0"><TR><TD>{1}<br/>{2}</TD></TR></TABLE>>, color="#000000", fillcolor="#ffffff", fontcolor="#000000"];\n'.format(node_idx, node_go_term, go_name)
			    else:
				line = lines[i]
			    
			    if node_go_term in go_color_pvalue_list:
				# GO term color
				color = go_color_pvalue_list[node_go_term]['color']
				line = line.replace(
				    'fillcolor="#ffffff"', 
				    'shape=rectangle, fillcolor="{0}"'.format(color))
				
				#line = line.replace('<br />', ' ')
				row = ''
				if 'pvalue' in go_color_pvalue_list[node_go_term]:
				    # List1 or List2 pvalue
				    p1 = go_color_pvalue_list[node_go_term]['pvalue']
				    row = '<TR><TD>P-value: {0}</TD></TR>'.format(p1)
				else:
				    # List1 only pvalue
				    p1 = go_color_pvalue_list[node_go_term]['pvalue_list1_o']
				    # Intersection pvalue
				    p2 = go_color_pvalue_list[node_go_term]['pvalue_inter']
				    # List2 only pvalue
				    p3 = go_color_pvalue_list[node_go_term]['pvalue_list2_o']
				    row = ''
				    
				    # List1 only
				    if labels[0] != '':
					row += '<TR><TD>P-value ({0}): {1}</TD></TR>'.format(labels[0].replace('_only', ''), p1) if p1 != '' else ''
				    
				    # Intersection
				    if labels[2] != '':
					row += '<TR><TD>P-value ({0}): {1}</TD></TR>'.format(labels[2].title(), p2) if p2 != '' else ''
				    
				    # List2 only
				    if labels[1] != '':
					row += '<TR><TD>P-value ({0}): {1}</TD></TR>'.format(labels[1].replace('_only', ''), p3) if p3 != '' else ''
				    
				line = line.replace('</TR>', '</TR>{0}'.format(row)).replace('<TABLE', '<TABLE WIDTH="150"')
			    else:
				line = line.replace(node_go_term+'<br/>', '').replace('fillcolor="#000000"', 'fillcolor="#ffffff"')
			    
			    line = line.replace(' CELLBORDER="1" CELLSPACING="0"', '')	
			else:
			    line = re.sub('color="#[0-9a-zA-Z]+"', 'color="#000000"', lines[i])			    
			    line = line.replace('arrowhead=none', 'arrowhead=normal')
			
			matchObj1 = re.match(r'^\s+node[0-9]+.+;$', line, re.M| re.I)
			if matchObj1:
			    o_f.write(line)
		    
		    o_f.write('}\n')
		    
		    # Add the i-th dot file
		    trees_list.append({
		        'dot_path': output_dot_file,
		        'pdf_path': output_pdf_file,
		        'work_dir_path': work_dir_path
		    })
		    
	return trees_list		    
		
		
	    