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
                 pathways_path, gene_2_pathway_path):
        
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
			            
        
    def executeEnrichmentAnalysis(self, prediction_dirs, no_cpu_cores, queue, 
                                  worker, resource_name):
        
        pred_dirs_per_core = [ [] for i in xrange(0, no_cpu_cores) ]
        
        # Prepare all prediction files to be assigned to each available core
        for i in xrange(0, len(prediction_dirs)):
            idx = int(math.fmod(i, no_cpu_cores))
            # Add the i-th prediction to the idx-th core
            pred_dirs_per_core[idx].append(prediction_dirs[i])
        
        jobs = []
        
        # For each core
        for i in xrange(0, len(pred_dirs_per_core)):
            if pred_dirs_per_core[i]:
                p = mp.Process(target = worker, 
                               args=(i,
                                     pred_dirs_per_core[i], # A pool of prediction(s) dir
                                     queue,
                                     self.refseq_gene_name_map, 
		                     self.pathways,
		                     self.gene_2_pathway,
		                     resource_name))               
                jobs.append(p)
                p.start() 
        
        # Wait worker(s) termination
        self.waitWorker(jobs, queue) 
        
                              
    def worker(self, idx, prediction_dirs, queue, refseq_gene_name_map, pathways,
               gene_2_pathway, resource_name):
        
        pred_json = {}
        pred_json_path = '/app/iso/static/public/output/prediction/prediction.json'
	
	if os.path.isfile(pred_json_path):
	    with open(pred_json_path, 'r') as f:
		pred_json = json.load(f)
				
        # For each miRNA predicted analysis
        for prediction in prediction_dirs:
            
            if pred_json:		
                predicted_targets_list = [] 
		dir_label = prediction['dir_label']
		mirna_id = None
		mir_pos = None
		mir_ref = None
		mir_mod = None
		mir_5p = None
		mir_3p = None
		
		matchObj1 = re.match(r'^(.+)__([0-9]+)_([a-z\+\-])_([a-z\+\-])$', 
	                            dir_label, 
	                            re.M|re.I)	
		
		matchObj2 = re.match(r'^(.+)__([0-9\+\-]+)_([0-9\+\-]+)$', 
	                            dir_label, 
	                            re.M|re.I)		
		if matchObj1:
		    # Substitution
		    mirna_id = matchObj1.group(1)
		    mir_pos = int(matchObj1.group(2))
		    mir_ref = matchObj1.group(3)
		    mir_mod = matchObj1.group(4)
		elif matchObj2:
		    # isomiR
		    mirna_id = matchObj2.group(1)
		    mir_5p = matchObj2.group(2)
		    mir_3p = matchObj2.group(3)	
		else:
		    # Wild Type
		    mirna_id = dir_label
		
		if mirna_id and mirna_id in pred_json:
		    preds = pred_json[mirna_id]
		    if mir_pos:
			# Substitution
			for s in preds['variations']['substitution']:
			    if s['position'] == mir_pos and \
			       s['ref'] == mir_ref and \
			       s['mod'] == mir_mod:
				predicted_targets_list = s['predictions']['targets']
				break
		    elif mir_5p:
			# isomiR
			for i in preds['variations']['isomir']:
			    if i['5p'] == mir_5p and \
			       i['3p'] == mir_3p:
				predicted_targets_list = i['predictions']['targets']
				break
		    else:
			predicted_targets_list = preds['predictions']['targets']
		
                # Enrichment 
                self.computeEnrichment(predicted_targets_list, 
                                       prediction['enrichment_dir_path'].replace('###', resource_name),
                                       refseq_gene_name_map, 
		                       pathways,
		                       gene_2_pathway)
		
	queue.put(({ 
	    'worker': idx
	}, idx))
            
            
    def comparisonWorker(self, idx, prediction_dirs, queue, refseq_gene_name_map, 
                         pathways, gene_2_pathway, resource_name):
	pred_json = {}
	pred_json_path = '/app/iso/static/public/output/prediction/prediction.json'
	
	if os.path.isfile(pred_json_path):
	    with open(pred_json_path, 'r') as f:
		pred_json = json.load(f)
		
        # For each miRNA predicted analysis
        for prediction in prediction_dirs:
            if 'wt' in prediction and \
               'variation' in prediction:
		
                # WT predictions dir label
                wt_pred_dir_label = prediction['wt']['dir_label']
                # Variation predictions dir label
                var_pred_dir_label = prediction['variation']['dir_label'] 
				
		# WT predicted targets list
		wt_pred_targets_list = [] 
		# Variation predicted targets list
		var_pred_targets_list = []     
		
		mirna_id = None
		mir_pos = None
		mir_ref = None
		mir_mod = None
		mir_5p = None
		mir_3p = None	
		
		matchObj1 = re.match(r'^(.+)__([0-9]+)_([a-z\+\-])_([a-z\+\-])$', 
	                            var_pred_dir_label, 
	                            re.M|re.I)	
		
		matchObj2 = re.match(r'^(.+)__([0-9\+\-]+)_([0-9\+\-]+)$', 
	                            var_pred_dir_label, 
	                            re.M|re.I)		
		if matchObj1:
		    # Substitution
		    mir_pos = int(matchObj1.group(2))
		    mir_ref = matchObj1.group(3)
		    mir_mod = matchObj1.group(4)
		elif matchObj2:
		    # isomiR
		    mir_5p = matchObj2.group(2)
		    mir_3p = matchObj2.group(3)
		
		if wt_pred_dir_label and wt_pred_dir_label in pred_json:
		    preds = pred_json[wt_pred_dir_label]
		    wt_pred_targets_list = preds['predictions']['targets']
		    if mir_pos:
			# Substitution
			for s in preds['variations']['substitution']:
			    if s['position'] == mir_pos and \
			       s['ref'] == mir_ref and \
			       s['mod'] == mir_mod:
				var_pred_targets_list = s['predictions']['targets']
				break
		    elif mir_5p:
			# isomiR
			for i in preds['variations']['isomir']:
			    if i['5p'] == mir_5p and \
			       i['3p'] == mir_3p:
				var_pred_targets_list = i['predictions']['targets']
				break
	
		wt_pred_targets_only = list(set(wt_pred_targets_list) - set(var_pred_targets_list))
		var_pred_targets_only = list(set(var_pred_targets_list) - set(wt_pred_targets_list))
		intersec_pred_targets = list(set(wt_pred_targets_list).intersection(set(var_pred_targets_list)))
		
		# Enrichment of the WT unique predicted targets
		self.computeEnrichment(wt_pred_targets_only, 
	                               prediction['wt_only_dir_path'].replace('###', resource_name),
	                               refseq_gene_name_map, 
		                       pathways,
		                       gene_2_pathway)
		
		# Enrichment of the variation unique predicted targets
		self.computeEnrichment(var_pred_targets_only, 
	                               prediction['variation_only_dir_path'].replace('###', resource_name),
	                               refseq_gene_name_map, 
		                       pathways,
		                       gene_2_pathway)
		
		# Enrichment of the WT and variation intersection
		self.computeEnrichment(intersec_pred_targets, 
	                               prediction['intersection_dir_path'].replace('###', resource_name),
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
		    
	    
    def computeEnrichment(self, predicted_targets_list, path, refseq_gene_name_map, 
                          pathways, gene_2_pathway):
	
	# If the directory exists
	if os.path.isdir(path):	
	    header = '\t'.join(['Pathway ID', 'Pathway name', 'P-value', 
		                '-Log(P-value)', 'Adjusted P-value', 
		                '-Log(Adjusted P-value)', 'Gene names']) + '\n'
	    
	    # ======================================================================
	    # TSV
	    # ======================================================================	
	    # Pathway enrichment file path
	    enrich_f_path = '%s/enrich.tsv' % path
	    enrich = []
	    
	    # ==================================================================
	    # JSON
	    # ==================================================================
	    # Pathway enrichment file path
	    enrich_f_json_path = '%s/enrich.json' % path
	    enrich_json = {}
	    
	    # Pathways list file path
	    enrich_f_pathway_json_path = '%s/enrich_pathway.json' % path
	    enrich_pathway_json = {}	
				    
	    if predicted_targets_list:
		input_genes = []    
		genes_pathways = []
		
		for refseq in predicted_targets_list:
		    if refseq in self.refseq_gene_name_map:
			gene_name = self.refseq_gene_name_map[refseq]['gene_name']
			
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
			
			# ======================================================
			# Fisher Exact Test
			# ======================================================
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
		    
		    # ==========================================================
		    # Benjamini/Hochberg (non-negative)
		    # ==========================================================
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
		    
		
