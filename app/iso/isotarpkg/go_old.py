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
        
    def __init__(self, nm_to_geneid_path, cmd_l, go_path):
        
        if not os.path.isfile(nm_to_geneid_path):
            self.ePrintClose("Error. No nm_to_geneid_path.tab file found.")
	
	if not os.path.isfile(go_path):
            self.ePrintClose("Error. No go.json file found.")

        self.debug = False          
        self.nm_to_geneid_path = nm_to_geneid_path
        self.refseq_gene_name_map = {}
        self.min_consensus = cmd_l.min_consensus
	self.pvalue = cmd_l.pvalue
	self.updated_dot_files_path_list = []
	self.go_path = go_path
	self.go_info = {}
	self.pvalue_type = cmd_l.pvalue_type
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
			    
			if gene_id not in self.geneid_to_genename:
			    self.geneid_to_genename[gene_id] = gene_name
        
        
    def executeEnrichmentAnalysis(self, prediction_dirs, no_cpu_cores, queue, 
                                  worker, obo_dag, geneid2gos_human, resource_name):
        
        self.getGeneIdToGeneNameList()
        
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
                                     obo_dag, 
		                     geneid2gos_human,
		                     resource_name))               
                jobs.append(p)
                p.start() 
        
        # Wait worker(s) termination
        self.waitWorker(jobs, queue) 
        
                              
    def worker(self, idx, prediction_dirs, queue, obo_dag, geneid2gos_human, 
               resource_name):
        
        pred_json = {}
        pred_json_path = '/app/iso/static/public/output/prediction/prediction.json'
	
	if os.path.isfile(pred_json_path):
	    with open(pred_json_path, 'r') as f:
		pred_json = json.load(f)
				
        # For each miRNA predicted analysis
        for prediction in prediction_dirs:
            
            if pred_json:
		filename_prefix = ''
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
		    
		    filename_prefix = '{0}_nt_var_{1}_{2}_{3}_'.format(
		        mirna_id, mir_pos, mir_ref, mir_mod)
		elif matchObj2:
		    # isomiR
		    mirna_id = matchObj2.group(1)
		    mir_5p = matchObj2.group(2)
		    mir_3p = matchObj2.group(3)	
		    filename_prefix = '{0}_isomir_{1}_{2}_'.format(
		        mirna_id, mir_5p, mir_3p)		    
		else:
		    # Wild Type
		    mirna_id = dir_label
		    filename_prefix = '{0}_wildtype_'.format(mirna_id)		    
		
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
                                       obo_dag, 
		                       geneid2gos_human)
		
		# GO terms Tree
		trees_list = self.run_r_go_term_dot_file(prediction['enrichment_dir_path'].replace(
		    '###', resource_name), filename_prefix)
		
		# Make the GO terms Tree
		for tree in trees_list:
		    self.create_go_terms_tree(tree['dot_path'], 
			                      tree['pdf_path'],
			                      tree['work_dir_path'])		
		
	queue.put(({ 
	    'worker': idx
	}, idx))
            
            
    def comparisonWorker(self, idx, prediction_dirs, queue, obo_dag, 
                         geneid2gos_human, resource_name):
	pred_json = {}
	pred_json_path = '/app/iso/static/public/output/prediction/prediction.json'
	
	if os.path.isfile(pred_json_path):
	    with open(pred_json_path, 'r') as f:
		pred_json = json.load(f)
		
        # For each miRNA predicted analysis
        for prediction in prediction_dirs:
            if 'wt' in prediction and \
               'variation' in prediction:
		
		filename_prefix = ''
		
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
		    mirna_id = matchObj1.group(1)
		    mir_pos = int(matchObj1.group(2))
		    mir_ref = matchObj1.group(3)
		    mir_mod = matchObj1.group(4)
		    
		    filename_prefix = '{0}_nt_var_{1}_{2}_{3}_comp_'.format(
		        mirna_id, mir_pos, mir_ref, mir_mod)		    
		elif matchObj2:
		    # isomiR
		    mirna_id = matchObj2.group(1)
		    mir_5p = matchObj2.group(2)
		    mir_3p = matchObj2.group(3)
		    filename_prefix = '{0}_isomir_{1}_{2}_comp_'.format(
		        mirna_id, mir_5p, mir_3p)		    
		else:
		    # Wild Type
		    mirna_id = wt_pred_dir_label
		
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
	                               obo_dag, 
	                               geneid2gos_human)
		
		# Enrichment of the variation unique predicted targets
		self.computeEnrichment(var_pred_targets_only, 
	                               prediction['variation_only_dir_path'].replace('###', resource_name),
	                               obo_dag, 
	                               geneid2gos_human)
		
		# Enrichment of the WT and variation intersection
		self.computeEnrichment(intersec_pred_targets, 
	                               prediction['intersection_dir_path'].replace('###', resource_name),
	                               obo_dag, 
	                               geneid2gos_human)   
		
		# GO terms Tree
		trees_list = self.run_r_go_term_comparison_tree(prediction['wt_only_dir_path'].replace('###', resource_name), 
	                                                        prediction['variation_only_dir_path'].replace('###', resource_name), 
	                                                        prediction['intersection_dir_path'].replace('###', resource_name),
		                                                filename_prefix)
		    
		# Make the GO terms Tree
		for tree in trees_list:
		    self.create_go_terms_tree(tree['dot_path'], 
			                      tree['pdf_path'],
			                      tree['work_dir_path'])		
		
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
		    
	    
    def computeEnrichment(self, predicted_targets_list, path, obo_dag, 
                          geneid2gos_human):
        
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
                            'Ratio in population', 'Gene names']) + '\n'
	
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
                
        if predicted_targets_list:
            input_genes = {}    
            
            for refseq in predicted_targets_list:
                if refseq in self.refseq_gene_name_map:
                    gene_name = self.refseq_gene_name_map[refseq]['gene_name']
                    gene_id = int(self.refseq_gene_name_map[refseq]['gene_id'])
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
				'go_namespace': 'cellular_component', 
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
	    if len(enrich_cc) > 0:
		f.write(''.join(enrich_cc))
	  
	with open(enrich_bp_f_path, 'w') as f:
	    f.write(header)
	    if len(enrich_bp) > 0:
		f.write(''.join(enrich_bp))
	  
	with open(enrich_mf_f_path, 'w') as f:
	    f.write(header)
	    if len(enrich_mf) > 0:
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
        
    
    
    def create_go_terms_dot_file(self, path, go_terms_list, work_dir_path):
	if len(go_terms_list) > 0:
	    command = '/usr/bin/Rscript /app/iso/isotarpkg/create_go_terms_tree.R --input {0} --go {1} --work {2}'.format(
	        path, ','.join(go_terms_list), work_dir_path)
	    i = 0
	    while i < self.give_me_my_dot_file:
		self.executeCommandLine(command) 
		
		if os.path.isfile(path):
		    break
		time.sleep(1)
		print 'Amigo 2 where are you?'
		i += 1
		
	    if not os.path.isfile(path):
		return False
	    
	    return True
	return False


    def create_go_terms_tree(self, dot_path, pdf_path, work_dir_path):
	if os.path.isfile(dot_path):
	    #time.sleep(2)
	    command = '/usr/bin/Rscript /app/iso/isotarpkg/create_updated_go_terms_tree.R --input {0} --output {1} --work {2}'.format(
	        dot_path, pdf_path, work_dir_path)
	    
	    self.executeCommandLine(command)  
	    
	    
    def run_r_go_term_dot_file(self, path, filename_prefix):
	
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
		    response = self.create_go_terms_dot_file(f_path, d.keys(), path)
		    
		    # GO terms color list
		    go_color_pvalue_list = { g: {
		        'color': '#1dc942',
		        'pvalue': d[g]['pvalue']
		    } for g in d.keys() }	
		    
		    if response:
			input_dot_file = '{0}/go_tree_{1}.dot'.format(path, go_c)
			output_dot_file = '{0}/go_tree_update_{1}.dot'.format(path, go_c)	
			output_pdf_file = '{0}/{1}go_tree_{2}.pdf'.format(path, 
			                                                  filename_prefix,
			                                                  go_c)
				
			# Update the Dot file and create te GO term Tree
			trees_list += self.update_go_tree_dot_file(input_dot_file, 
			                                           output_dot_file,
			                                           output_pdf_file,
			                                           go_color_pvalue_list,
			                                           path)
		    
	return trees_list
		    
		    
    def run_r_go_term_comparison_tree(self, wt_o_path, var_o_path, inter_path, 
                                      filename_prefix):
	
	# Comparison directory path
	comp_path = wt_o_path.replace('/wt_only', '')	
	go_class = ['bp', 'cc', 'mf']
	trees_list = []
		
	colors = {
	    '001': '#ff6347', # Variation: Red
	    '010': '#3cb371', # Intersection: Green
	    '011': '#9e8b5c', # Intersection & Variation:Brown
	    '100': '#1e90ff', # Wild Type: Blue
	    '101': '#8f7aa3', # Wild Type & Variation: Violet
	    '110': '#2da2b8', # Wild Type & Intersection: Dodger Blue
	    '111': '#b4b4b4' # All: Gray
	}
	    
	# Extract from each JSON the GO terms list
        for go_c in go_class:
	    go_list = []
	    comp = {
		'wt_only': [], 
		'intersection': [], 
		'variation_only': []
	    }	
	    
	    comp_pvalue = {
		'wt_only': {}, 
		'intersection': {}, 
		'variation_only': {}
	    }	    
	    
	    # Retrieve GO terms from Wild Type only, Intersection, and Variation only
	    for c, l in comp.iteritems():
		enrich_f_json_path = '{0}/{1}/enrich_{2}.json'.format(comp_path, 
		                                                      c, 
		                                                      go_c)
		
		if os.path.isfile(enrich_f_json_path):
		    with open(enrich_f_json_path, 'r') as f:
			d = json.load(f)
			comp[c] = d.keys()			
			comp_pvalue[c] = { g: d[g]['pvalue'] for g in d.keys() }
	    
	    # Unique GO term list
	    go_list = list(set(comp['wt_only']+comp['intersection']+comp['variation_only']))
	    
	    # Dot file path
	    f_path = '{0}/go_tree_{1}.dot'.format(comp_path, go_c)
	    
	    # Make the GO terms Tree
	    response = self.create_go_terms_dot_file(f_path, go_list, comp_path)
	    
	    go_color_pvalue_list = {}
	    
	    # For each GO term checks in which list the GO  term is
	    for g in go_list:
		b1 = '1' if g in comp['wt_only'] else '0'
		b2 = '1' if g in comp['intersection'] else '0'
		b3 = '1' if g in comp['variation_only'] else '0'
		
		label = '{0}{1}{2}'.format(b1, b2, b3)
		   
		if label in colors:
		    # GO terms color list
		    go_color_pvalue_list[g] = {
		        'color': colors[label],
		        'pvalue_wt_o': comp_pvalue['wt_only'][g] if g in comp_pvalue['wt_only'] else '',
		        'pvalue_inter': comp_pvalue['intersection'][g] if g in comp_pvalue['intersection'] else '',
		        'pvalue_var_o': comp_pvalue['variation_only'][g] if g in comp_pvalue['variation_only'] else ''
		    }
	    
	    if response:
		input_dot_file = '{0}/go_tree_{1}.dot'.format(comp_path, go_c)
		output_dot_file = '{0}/go_tree_update_{1}.dot'.format(comp_path, go_c)	
		output_pdf_file = '{0}/{1}go_tree_{2}.pdf'.format(comp_path, 
		                                                  filename_prefix,
		                                                  go_c)
				
		# Update the Dot file and create te GO term Tree
		trees_list += self.update_go_tree_dot_file(input_dot_file, 
		                                           output_dot_file,
		                                           output_pdf_file,
		                                           go_color_pvalue_list,
		                                           comp_path)
	return trees_list
	    
    
    def update_go_tree_dot_file(self, input_dot_file, output_dot_file,
                                output_pdf_file, go_color_pvalue_list,
                                work_dir_path):
	trees_list = []
	
	merged_lines = []
	
	if os.path.isfile(input_dot_file) and len(go_color_pvalue_list) > 0:
	    with open(output_dot_file, 'w') as o_f:
		with open(input_dot_file, 'r') as i_f:
		    lines = i_f.readlines()
		    
		    # //////////////////////////////////////////////////////////
		    # Merge all splitted lines for each node/edge
		    # //////////////////////////////////////////////////////////
		    '''
		    row = ''
		    for i in xrange(1, len(lines)):
			matchObj = re.match(r',$', lines[i], re.M| re.I)
			if matchObj:
			    row += lines[i].strip()
			else:
			    if row == '':
				row = lines[i].strip()
			    merged_lines.append(row + '\n')
			    row = ''
		    
		    lines = merged_lines
		    '''
		    
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
				    while not re.match(r'^\s+node[0-9]+\s+\[.+TD\>([a-z]+:[a-z0-9]+).+\];$', lines[i], re.M| re.I):
					i += 1
					
				    # For the next node
				    i -= 1
					
				if node_go_term in self.go_info:
				    go_name = self.go_info[node_go_term]['name']
				    # Replace bad characters
				    go_name = go_name.replace('&', '&amp;').replace('"', '&#34;').replace('<', '&lt;;').replace('>', '&gt;')
				    
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
				    # Wild Type or Variation pvalue
				    p1 = go_color_pvalue_list[node_go_term]['pvalue']
				    row = '<TR><TD>P-value: {0}</TD></TR>'.format(p1)
				else:
				    # Wild Type only pvalue
				    p1 = go_color_pvalue_list[node_go_term]['pvalue_wt_o']
				    # Intersection pvalue
				    p2 = go_color_pvalue_list[node_go_term]['pvalue_inter']
				    # Variation only pvalue
				    p3 = go_color_pvalue_list[node_go_term]['pvalue_var_o']
				    row = '<TR><TD>P-value (Wild Type): {0}</TD></TR>'.format(p1) if p1 != '' else ''
				    row += '<TR><TD>P-value (Intersection): {0}</TD></TR>'.format(p2) if p2 != '' else ''
				    row += '<TR><TD>P-value (Variation): {0}</TD></TR>'.format(p3) if p3 != '' else ''
				    
				line = line.replace('</TR>', '</TR>{0}'.format(row)).replace('<TABLE', '<TABLE WIDTH="130"')
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
		    	
			
    #/usr/bin/Rscript /app/iso/isotarpkg/create_go_terms_tree.R --path /output/enrichment/hsa-let-7c-5p__hsa-let-7c/nt_sub/hsa-let-7c-5p__17_A_G --class "bp" --stage "tree_file"
		    
		
		
	    