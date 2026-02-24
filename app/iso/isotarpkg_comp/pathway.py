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
from scipy import stats
from statsmodels.stats.multitest import multipletests

class Pathway():
        
    def __init__(self, nm_to_geneid_path, pathway_resource_name, cmd_l, 
                 pathways_path, gene_2_pathway_path, output_path):
        
        # Supported pathway DBs
        self.pathway_reources_list = ['kegg', 'reactome']
	
        if pathway_resource_name not in self.pathway_reources_list:
	    self.ePrintClose("Error. A wrong pathway DB has been specified. Please choose among: {0}".format(
	        ', '.join(self.pathway_reources_list)))
        
        if not os.path.isfile(pathways_path):
	    self.ePrintClose("Error. No pathway.json file found.")
	
	if not os.path.isfile(gene_2_pathway_path):
	    self.ePrintClose("Error. No genename_pathways.json file found.")
	    
	if not os.path.isfile(nm_to_geneid_path):
            self.ePrintClose("Error. No nm_to_geneid_path.tab file found.")
	    
	if not os.path.isdir(output_path):
            self.ePrintClose("Error. No enrichment_comp directory found.")
	    
	self.output_enrichment_dir = output_path
	self.pathways_path = pathways_path
	self.pvalue = cmd_l.pvalue
	self.pvalue_type = cmd_l.pvalue_type
        self.gene_2_pathway_path = gene_2_pathway_path
	self.nm_to_geneid_path = nm_to_geneid_path
	self.resource_name = pathway_resource_name

        self.debug = False          
        self.gene_2_pathway = {}
	self.pathways = {}
	self.refseq_gene_name_map = {}
	self.genes_list = {}
	self.genes_univ_list = []
	self.input_targets_list1 = [ el.strip() for el in cmd_l.targets_list1.replace('"', '').replace('\'', '').split(',') if el.strip() != '' ]  
	self.input_targets_list2 = [ el.strip() for el in cmd_l.targets_list2.replace('"', '').replace('\'', '').split(',') if el.strip() != '' ]	
	self.discarded_genes_in_list1 = []
	self.discarded_genes_in_list2 = []
	self.targets_list1 = []  
	self.targets_list2 = []
	self.name_list1 = re.sub(r'\s+', '_', re.sub(r'[^0-9a-zA-Z_\s\-]+', '', cmd_l.name_list1))  
	self.name_list2 = re.sub(r'\s+', '_', re.sub(r'[^0-9a-zA-Z_\s\-]+', '', cmd_l.name_list2)) 	
	self.support = {}
	
	with open(self.gene_2_pathway_path, 'r') as f:
	    self.gene_2_pathway = json.load(f)
	    
	with open(self.pathways_path, 'r') as f:
	    self.pathways = json.load(f)
	    
	self.getGeneIdToGeneNameList()
	
	
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
    
                  
    def getGeneIdToGeneNameList(self): 
        
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
			
			if gene_name not in self.genes_list:			    
                            self.genes_list[gene_name] = {
                                'gene_name': gene_name,
                                'pathways': self.gene_2_pathway[self.resource_name][gene_name]['pathways'] 
			        if gene_name in self.gene_2_pathway[self.resource_name] else []	
                            }
			
			
        self.genes_univ_list = self.genes_list.keys()
	
    
    def checkInputGenesList(self): 
		
	# Check the first genes list
        for gene in self.input_targets_list1:
	    gene_name = ''
	    if gene in self.refseq_gene_name_map:
		g_name = self.refseq_gene_name_map[gene]['gene_name']
		if g_name in self.genes_list:
		    gene_name = g_name
	    elif gene in self.genes_list:
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
		if g_name in self.genes_list:
		    gene_name = g_name
	    elif gene in self.genes_list:
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
			            
        
    def executeEnrichmentAnalysis(self, no_cpu_cores, queue, worker, 
                                  resource_name, utils):
	
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
                p = mp.Process(target = worker, 
                               args=(i,
                                     enrich_dirs_per_core[i], # A pool of prediction(s) dir
                                     queue,
                                     self.refseq_gene_name_map, 
		                     self.pathways,
		                     self.gene_2_pathway,
		                     resource_name))               
                jobs.append(p)
                p.start() 
        
        # Wait worker(s) termination
        self.waitWorker(jobs, queue) 
        
                              
    def worker(self, idx, enrichment_dirs, queue, refseq_gene_name_map, pathways,
               gene_2_pathway, resource_name):
        
        # For each enrichment dir
        for enrichment in enrichment_dirs:
	    
            if len(enrichment['targets']) > 0:
                # Enrichment
                self.computeEnrichment(enrichment['targets'], 
                                       enrichment['enrichment_dir_path'],
                                       refseq_gene_name_map, 
		                       pathways,
		                       gene_2_pathway)
		
	queue.put(({ 
	    'worker': idx
	}, idx))
                 
    
    def waitWorker(self, jobs, queue):
        
        # Collects all results into in a list                    
        for j in xrange(0, len(jobs)):
            # Get worker result
            res = queue.get()
            if res:
                worker_id = res[0]['worker']
                
        # Wait for all worker processes to finish
        for p in jobs:
            p.join() 
		    
	    
    def computeEnrichment(self, targets_list, path, refseq_gene_name_map, 
                          pathways, gene_2_pathway):
		
        header = '\t'.join(['Pathway ID', 'Pathway name', 'P-value', 
                            '-Log(P-value)', 'Adjusted P-value', 
                            '-Log(Adjusted P-value)', 'Gene names']) + '\n'
	
	# ======================================================================
	# TSV
	# ======================================================================	
        # Pathway enrichment file path
        enrich_f_path = '%s/enrich.tsv' % path
	enrich = []
	
	# ======================================================================
	# JSON
	# ======================================================================
	# Pathway enrichment file path
	enrich_f_json_path = '%s/enrich.json' % path
	enrich_json = {}
	
	# Pathways list file path
	enrich_f_pathway_json_path = '%s/enrich_pathway.json' % path
	enrich_pathway_json = {}	
		                
        if targets_list:
            input_genes = []    
            genes_pathways = [] 
	    
	    for gene_name in targets_list:
		# Prepare the gene names list
		if gene_name not in input_genes:
		    input_genes.append(gene_name)
		
		# Get all the pathways associated with each gene name into input_genes list
		if gene_name in self.genes_list:
		    genes_pathways += [ pathway_id for pathway_id in 
		                        self.genes_list[gene_name]['pathways'] 
		                        if pathway_id not in genes_pathways ]
	           
	    # For each pathway ID into genes_pathways list
	    if genes_pathways:
		pvalues_list = []
		p_tmp = []
		
		for pathway_id in genes_pathways:
		    # All available gene names
		    m1_1 = len(list(set(input_genes).intersection(self.pathways[self.resource_name][pathway_id]['genes'])))
		    m1_2 = len(list(set(self.pathways[self.resource_name][pathway_id]['genes']) - set(input_genes)))
		    m2_1 = len(list(set(input_genes) - set(self.pathways[self.resource_name][pathway_id]['genes'])))
		    m2_2 = len(list(set(self.genes_univ_list) - set(input_genes).union(self.pathways[self.resource_name][pathway_id]['genes'])))
		    
		    obs = [[m1_1, m1_2],[m2_1, m2_2]]
		    		    
		    # ==========================================================
		    # Fisher Exact Test
		    # ==========================================================
		    oddsratio, p_uncorrected = stats.fisher_exact(obs)
		    pvalues_list.append(p_uncorrected)
		    
		    l = [
			pathway_id, 
			self.pathways[self.resource_name][pathway_id]['pathway_name'], 
			'{0:5.2e}'.format(p_uncorrected), 
			'{:f}'.format(math.fabs(math.log10(p_uncorrected))), 
			'', 
			'',
			'{0}'.format(';'.join(list(set(input_genes).intersection(self.pathways[self.resource_name][pathway_id]['genes']))))
		    ]	
		    
		    p_tmp.append(l)		    
		
		# ==============================================================
		# Benjamini/Hochberg (non-negative)
		# ==============================================================
		reject, pvals_corrected, alphacSidak, alphacBonf = multipletests(pvalues_list, method='fdr_bh')
		
		# Add all adjusted pvalue
		for i in xrange(0, len(p_tmp)):
		    if i < len(pvals_corrected):
			p_tmp[i][4] = '{0:5.2e}'.format(pvals_corrected[i])
			p_tmp[i][5] = '{:f}'.format(math.fabs(math.log10(pvals_corrected[i])))
			
			enrich.append('\t'.join(p_tmp[i]).encode('utf-8') + '\n') 
			enrich_pathway_json[p_tmp[i][0]] = math.fabs(math.log10(float(p_tmp[i][2])))	
			
			# Use a specific pvalue type: raw or adj
			pvalue = float(p_tmp[i][2])
			if self.pvalue_type == 'adj':
			    pvalue = float(pvals_corrected[i])
			    
			# P-value filtering according to the user
			if pvalue < self.pvalue:
			    d = {
				'pathway_name': p_tmp[i][1], 
				'pathway_id': p_tmp[i][0], 
				'pvalue': p_tmp[i][2], 
				'-Log(pvalue)': p_tmp[i][3], 
				'adj_pvalue': p_tmp[i][4], 
				'-Log(adj_pvalue)': p_tmp[i][5],
				'targets_list': p_tmp[i][6]
			    }
			    
			    enrich_json[p_tmp[i][0]] = d
			    
	# ==================================================================
        # TSV
	# ==================================================================
	with open(enrich_f_path, 'w') as f:
	    f.write(header)
	    if len(enrich) > 0:
		f.write(''.join(enrich))
	    
	# ==================================================================
	# JSON
	# ==================================================================
	with open(enrich_f_json_path, 'w') as f:
	    f.write(json.dumps(enrich_json, indent=4)) 
	    
	# ==================================================================
	# JSON Pathway list
	# ==================================================================
	with open(enrich_f_pathway_json_path, 'w') as f:
	    f.write(json.dumps(enrich_pathway_json, indent=4)) 	
		    
		
