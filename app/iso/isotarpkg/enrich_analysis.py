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
import go as go
import pathway as pt
import json
import csv
import numpy as np
from goatools.base import download_go_basic_obo
from goatools.base import download_ncbi_associations
from goatools.obo_parser import GODag
from goatools.associations import read_ncbi_gene2go
from goatools.test_data.genes_NCBI_9606_ProteinCoding import GeneID2nt as GeneID2nt_hsa
from goatools.go_enrichment import GOEnrichmentStudy
from operator import itemgetter
import math

class EnrichmentAnalysis():
        
    def __init__(self, output_dir, nm_to_geneid_path, cmd_l, obo_dag_path, 
                 gene_2_go_path, go_path, pathways_path, gene_2_pathway_path):
        
        self.output_dir = output_dir
        self.output_prediction_dir = '{0}/prediction'.format(output_dir)
        self.output_enrichment_dir = '{0}/enrichment'.format(output_dir)
        self.nm_to_geneid_path = nm_to_geneid_path
	self.obo_dag_path = obo_dag_path
	self.gene_2_go_path = gene_2_go_path
	self.pvalue = cmd_l.pvalue
	self.go_path = go_path
	self.pathways_path = pathways_path
	self.gene_2_pathway_path = gene_2_pathway_path	
	self.enrich_type = cmd_l.enrich_type
            
        if not os.path.isdir(self.output_dir):
            self.ePrintClose("Error. No output directory found.")
        
        if not os.path.isdir(self.output_prediction_dir):
            self.ePrintClose("Error. No output prediction directory found.")
        
        if not os.path.isfile(self.nm_to_geneid_path):
            self.ePrintClose("Error. No nm_to_geneid_path.tab file found.")
	
	if not os.path.isfile(self.obo_dag_path):
	    self.ePrintClose("Error. No go-basic.obo file found.")
	
	if not os.path.isfile(self.gene_2_go_path):
	    self.ePrintClose("Error. No gene2go file found.")
	    
	if not os.path.isfile(self.go_path):
            self.ePrintClose("Error. No go.json file found.")
	    
	if not os.path.isfile(self.pathways_path):
	    self.ePrintClose("Error. No pathways.json file found.")
	
	if not os.path.isfile(self.gene_2_pathway_path):
	    self.ePrintClose("Error. No genename_pathways.json file found.")
            
        # The object provides a set of facilities to create and remove the output
        # directory, as well as functions to run bash commands and multi processes
        utl = ut.ITUtils(cmd_l.no_cpu_cores)
        
        self.enrich_dirs_list = []
        
        # Used to store output.tsv
        self.enrich_comparison_dirs_list = []
        
        # Get all miRNA(s) prediction directories
        predictions = next(os.walk(self.output_prediction_dir))[1] 
        
        # Delete the enrichment directory
        utl.removeDirectory(self.output_enrichment_dir)
	utl.createDirectory(self.output_enrichment_dir)
                        
        for prediction in predictions:
            mirna_json_path = '{0}/{1}/mirna_check.json'.format(
	        self.output_prediction_dir, 
                prediction) 
	    
	    # miRNA root direcotry
	    mirna_root_path = '{0}/{1}'.format(self.output_prediction_dir, 
	                                       prediction) 
	    
            if os.path.isfile(mirna_json_path):
                with open(mirna_json_path, 'r') as f:
                    mirna_info = json.load(f)
                    
                    if mirna_info and \
                       'identifier' in mirna_info and \
                       'output_path' in mirna_info:    
                        enrich_root_dir = self.prepareEnrichedFiles(
                            mirna_info['identifier'], 
                            mirna_info['output_path'],
			    prediction)
                        
                        if enrich_root_dir:
                            self.enrich_dirs_list.append(enrich_root_dir) 
			    
			    # Create all root enrichment type directories
			    for en_t in self.enrich_type:			    
				# Create a new empty output directory
				utl.createDirectory(enrich_root_dir['enrichment_dir_path'].replace('###', en_t))
				
			    available_variation = False
                
                            # ==================================================
                            # ==================================================
                            # Variations 
                            # ==================================================
                            # ==================================================
                            if 'variations' in mirna_info:
                                if 'substitution' in mirna_info['variations']:				    
                                    # Nucleotide substitution
                                    for j in xrange(0, len(
                                        mirna_info['variations']['substitution'])):
					available_variation = True
                                        mir = mirna_info['variations']['substitution'][j]
                                        
                                        enrich_dir = self.prepareEnrichedFiles(
                                            mir['identifier'], 
                                            mir['output_path'],
					    prediction)
					
                                        if enrich_dir:
                                            self.enrich_dirs_list.append(
                                                enrich_dir) 
					    
					    # Create all root enrichment type directories
					    for enrich_opt in self.enrich_type:			    
						# Create a new empty output directory
						utl.createDirectory(enrich_dir['enrichment_dir_path'].replace('###', enrich_opt))
						
                                            # ==================================
                                            # WT vs. variation 
                                            # ==================================                                     
                                            p = '{0}/wt_vs_variation'.format(
					        enrich_dir['enrichment_dir_path'])
                                            
                                            p1 = '{0}/wt_only'.format(p)
                                            p2 = '{0}/intersection'.format(p)
                                            p3 = '{0}/variation_only'.format(p)
                                            
					    # Create all enrichment type directories
					    for enrich_opt in self.enrich_type:			    						
						utl.createDirectory(p1.replace('###', enrich_opt))                                        
						utl.createDirectory(p2.replace('###', enrich_opt))
						utl.createDirectory(p3.replace('###', enrich_opt))
                                            
                                            self.enrich_comparison_dirs_list.append({
                                                'wt': {
                                                    'predictions_file_path': 
                                                    enrich_root_dir['predictions_file_path'],
                                                    'enrichment_dir_path': 
                                                    enrich_root_dir['enrichment_dir_path'],
                                                    'dir_label': 
                                                    enrich_root_dir['dir_label']
                                                },
                                                'variation': {
                                                    'predictions_file_path': 
                                                    enrich_dir['predictions_file_path'],
                                                    'enrichment_dir_path': 
                                                    enrich_dir['enrichment_dir_path'],
                                                    'dir_label': 
                                                    enrich_dir['dir_label']
                                                },
                                                'wt_only_dir_path': p1,
                                                'intersection_dir_path': p2,
                                                'variation_only_dir_path': p3
                                            })
                            
                            if 'variations' in mirna_info:
                                if 'isomir' in mirna_info['variations']:				    
                                    # isomiR      
                                    for j in xrange(0, len(
                                        mirna_info['variations']['isomir'])):
					available_variation = True
                                        mir = mirna_info['variations']['isomir'][j]
                                        
                                        enrich_dir = self.prepareEnrichedFiles(
                                            mir['identifier'], 
                                            mir['output_path'],
					    prediction)
					
                                        if enrich_dir:
                                            self.enrich_dirs_list.append(
                                                enrich_dir) 
                                            
					    # Create all root enrichment type directories
					    for enrich_opt in self.enrich_type:			    
						# Create a new empty output directory
						utl.createDirectory(enrich_dir['enrichment_dir_path'].replace('###', enrich_opt))                                            
                                            
                                            # ==================================
                                            # WT vs. variation 
                                            # ==================================                                     
                                            p = '{0}/wt_vs_variation'.format(
					        enrich_dir['enrichment_dir_path'])

                                            p1 = '{0}/wt_only'.format(p)
                                            p2 = '{0}/intersection'.format(p)
                                            p3 = '{0}/variation_only'.format(p)
                                            
					    # Create all enrichment type directories
					    for enrich_opt in self.enrich_type:			    						
						utl.createDirectory(p1.replace('###', enrich_opt))                                        
						utl.createDirectory(p2.replace('###', enrich_opt))
						utl.createDirectory(p3.replace('###', enrich_opt))  
                                            					    
                                            self.enrich_comparison_dirs_list.append({
                                                'wt': {
                                                    'predictions_file_path': 
                                                    enrich_root_dir['predictions_file_path'],
                                                    'enrichment_dir_path': 
                                                    enrich_root_dir['enrichment_dir_path'],
                                                    'dir_label': 
                                                    enrich_root_dir['dir_label']
                                                },
                                                'variation': {
                                                    'predictions_file_path': 
                                                    enrich_dir['predictions_file_path'],
                                                    'enrichment_dir_path': 
                                                    enrich_dir['enrichment_dir_path'],
                                                    'dir_label': 
                                                    enrich_dir['dir_label']
                                                },
                                                'wt_only_dir_path': p1,
                                                'intersection_dir_path': p2,
                                                'variation_only_dir_path': p3
                                            })  
			    
			    # Only for those miRNA without any variation
                            if not available_variation:
				self.enrich_comparison_dirs_list.append({
				    'wt': {
				        'predictions_file_path': 
				        enrich_root_dir['predictions_file_path'],
				        'enrichment_dir_path': 
				        enrich_root_dir['enrichment_dir_path'],
				        'dir_label': 
				        enrich_root_dir['dir_label']
				    },
				    'variation': {
				        'predictions_file_path': '',
				        'enrichment_dir_path': '',
				        'dir_label': ''
				    },
				    'wt_only_dir_path': '',
				    'intersection_dir_path': '',
				    'variation_only_dir_path': ''
				})
        
        if self.enrich_dirs_list:   
	    
	    for enrich_opt in self.enrich_type:
		
		# ==============================================================
		# GO term Enrichment Analysis
		# ==============================================================
		if 'go' == enrich_opt:
		    obo_dag = GODag(self.obo_dag_path)
		    geneid2gos_human = read_ncbi_gene2go(self.gene_2_go_path, 
			                                 taxids=[9606])
		    # Workers queue
		    go_queue = mp.Queue()
	    
		    # Set the GO term object
		    GO = go.GO(self.nm_to_geneid_path, cmd_l, self.go_path)
				
		    trees_list = []
		    trees_l1 = []
		    trees_l2 = []
		    		    
		    # Execute the go term enrichment analysis for all miRNA(s) 
		    GO.executeEnrichmentAnalysis(self.enrich_dirs_list, 
			                         cmd_l.no_cpu_cores, 
			                         go_queue,
			                         GO.worker,
			                         obo_dag, 
			                         geneid2gos_human,
		                                 enrich_opt)
		    # Workers queue
		    go_queue = mp.Queue()
		    
		    if self.enrich_comparison_dirs_list:
			# Execute the wt vs. variation go term enrichment analysis for all miRNA(s) 
			GO.executeEnrichmentAnalysis(self.enrich_comparison_dirs_list, 
				                     cmd_l.no_cpu_cores, 
				                     go_queue,
				                     GO.comparisonWorker,
				                     obo_dag, 
				                     geneid2gos_human,
			                             enrich_opt)
		
		else:
		    # ==========================================================
		    # Pathway Enrichment Analysis
		    # ==========================================================		    
		    # Set the Pathway object
		    Pathway = pt.Pathway(self.nm_to_geneid_path, enrich_opt, 
		                         cmd_l, self.pathways_path,
		                         self.gene_2_pathway_path)
		    
		    # Workers queue
		    pt_queue = mp.Queue()
		    
		    # Execute the pathway enrichment analysis for all miRNA(s) 
		    Pathway.executeEnrichmentAnalysis(self.enrich_dirs_list, 
		                                      cmd_l.no_cpu_cores, 
		                                      pt_queue,
		                                      Pathway.worker,
		                                      enrich_opt)
		    # Workers queue
		    pt_queue = mp.Queue()
		    
		    if self.enrich_comparison_dirs_list:
			# Execute the pathway enrichment analysis for all miRNA(s) 
			Pathway.executeEnrichmentAnalysis(self.enrich_comparison_dirs_list, 
			                                  cmd_l.no_cpu_cores, 
			                                  pt_queue,
			                                  Pathway.comparisonWorker,
			                                  enrich_opt)	
			
		self.prepareEnrichedJsonFile(self.enrich_comparison_dirs_list)
	             
                
    def prepareEnrichedFiles(self, mirna, path, dirname):
        
        enrich_dir = {}
        predictions_file_path = '{0}/output.tsv'.format(path)
        
        if os.path.isfile(predictions_file_path):
            e_path = path.replace('prediction', 'enrichment').replace(dirname, '{0}/###'.format(dirname))
            enrich_dir = {
               'dir_label': mirna,
               'enrichment_dir_path': e_path,
               'predictions_file_path': predictions_file_path
            }
            
        return enrich_dir
        
        
    def prepareEnrichedJsonFile(self, enrich_dirs_list):
        
	en_path = '{0}/enrichment.json'.format(self.output_enrichment_dir)
	master_enrich_json = {}
	
        for enrich_opt in self.enrich_type:
	    enrich_json = {}
	    
	    # GO term
	    if enrich_opt == 'go':        
		enrich_files = {
		    'bp': 'enrich_bp.json',
		    'cc': 'enrich_cc.json',
		    'mf': 'enrich_mf.json'
		}
		
		# GO -Log(pvalues)
		go_pvalues_files = {
		    'bp': 'enrich_bp_go.json',
		    'cc': 'enrich_cc_go.json',
		    'mf': 'enrich_mf_go.json'
		}	
			
		for mir in enrich_dirs_list: 
		    filename_prefix = ''
		    mirna_id = mir['wt']['dir_label']
		    wt_path = mir['wt']['enrichment_dir_path'].replace('###', enrich_opt)
		    
		    var_id = mir['variation']['dir_label'].split('__')
		    var_path = mir['variation']['enrichment_dir_path'].replace('###', enrich_opt)    
		    
		    # ==================================================================
		    # ==================================================================
		    # Get all Wild Type GO term pvalues
		    # ==================================================================
		    # ==================================================================
		    wt_go_pvalue = { 'bp': [], 'cc': [], 'mf': [] }
		    var_go_pvalue = { 'bp': [], 'cc': [], 'mf': [] }
		    histo_go_pvalue = { 'bp': [], 'cc': [], 'mf': [] }	    
		    
		    for go_class, go_file in go_pvalues_files.iteritems():
			f_path = '{0}/{1}'.format(wt_path, go_file)
			if os.path.isfile(f_path):
			    with open(f_path, 'r') as f:
				d = json.load(f)
				wt_go_pvalue[go_class] = d	
		    
		    for go_class, go_file in go_pvalues_files.iteritems():
			f_path = '{0}/{1}'.format(var_path, go_file)
			if os.path.isfile(f_path):
			    with open(f_path, 'r') as f:
				d = json.load(f)
				var_go_pvalue[go_class] = d
		    
		    # //////////////////////////////////////////////////////////////////
		    # Wild Type
		    # //////////////////////////////////////////////////////////////////
		    if mirna_id not in enrich_json:
			enrich_json[mirna_id] = {
			    'mirna_id': mirna_id,
			    'variations': {
				'substitution': [],
				'isomir': []
			    },
			    'bp': {},
			    'cc': {},
			    'mf': {},
			    'threshold': math.fabs(math.log10(self.pvalue)),
			    'bp_tree': '',
			    'cc_tree': '',
			    'mf_tree': ''
			}
			
		    
		    # Get all Wild Type GO terms for each GO class
		    for go_class, go_file in enrich_files.iteritems():
			f_path = '{0}/{1}'.format(wt_path, go_file) 
			if os.path.isfile(f_path):
			    with open(f_path, 'r') as f:
				d = json.load(f)
				enrich_json[mirna_id][go_class] = d
				
				filename_prefix = '{0}_wildtype_'.format(mirna_id)				
				
				# Wild Type GO Tree
				tree_path = '{0}/{1}go_tree_{2}.pdf'.format(wt_path, 
				                                            filename_prefix,
				                                            go_class) 
				tree_label = '{0}_tree'.format(go_class)
				if os.path.isfile(tree_path):
				    enrich_json[mirna_id][tree_label] = 'public{0}'.format(tree_path.replace('/app/iso/static/public', ''))
									   
		    # //////////////////////////////////////////////////////////////////
		    # Get all variation GO terms for each GO class               
		    # //////////////////////////////////////////////////////////////////
		    if mir['variation']['dir_label'] != '':
			variation = {}
			
			var_type = None
			matchObj1 = re.match(r'^([0-9]+)_([^_]+)_([^_]+)$', var_id[1], re.M|re.I)
			matchObj2 = re.match(r'^([^_]+)_([^_]+)$', var_id[1], re.M|re.I)
			
			if matchObj1:  
			    var_type = 'substitution'
			    variation['position'] = int(matchObj1.group(1))
			    variation['ref'] = matchObj1.group(2)
			    variation['mod'] = matchObj1.group(3)
			    filename_prefix = '{0}_nt_var_{1}_{2}_{3}_'.format(
				mirna_id, 
			        variation['position'], 
			        variation['ref'], 
			        variation['mod'])			    
			elif matchObj2:
			    var_type = 'isomir'
			    variation['5p'] = matchObj2.group(1)
			    variation['3p'] = matchObj2.group(2)  
			    filename_prefix = '{0}_isomir_{1}_{2}_'.format(
				mirna_id, variation['5p'], variation['3p'])			    
		
			if var_type:
			    variation['bp'] = {}
			    variation['cc'] = {}
			    variation['mf'] = {}
			    variation['wt_only'] = {}
			    variation['intersection'] = {}
			    variation['variation_only'] = {}  
			    variation['bp_his'] = []
			    variation['cc_his'] = []
			    variation['mf_his'] = []
			    variation['comp_bp_his'] = []
			    variation['comp_cc_his'] = []
			    variation['comp_mf_his'] = []
			    variation['bp_tree'] = ''
			    variation['cc_tree'] = ''
			    variation['mf_tree'] = ''
			    variation['bp_wt_vs_variation_tree'] = ''
			    variation['cc_wt_vs_variation_tree'] = ''
			    variation['mf_wt_vs_variation_tree'] = ''		
			    
			    # Get all Variation GO terms for each GO class
			    for go_class, go_file in enrich_files.iteritems():
				f_path = '{0}/{1}'.format(var_path, go_file) 
				if os.path.isfile(f_path):
				    with open(f_path, 'r') as f:
					d = json.load(f)
					variation[go_class] = d  
					wt_go = enrich_json[mirna_id][go_class].keys()
					histo_go_pvalue[go_class] = list(set(wt_go + d.keys()))
					
					# Variation GO Tree
					tree_path = '{0}/{1}go_tree_{2}.pdf'.format(
					    var_path, filename_prefix, go_class) 
					tree_label = '{0}_tree'.format(go_class)
					if os.path.isfile(tree_path):
					    variation[tree_label] = 'public{0}'.format(tree_path.replace('/app/iso/static/public', ''))
					    
					# Comparison GO Tree
					tree_path = '{0}/wt_vs_variation/{1}comp_go_tree_{2}.pdf'.format(
					    var_path, filename_prefix, go_class) 
					tree_label = '{0}_wt_vs_variation_tree'.format(go_class)
					if os.path.isfile(tree_path):
					    variation[tree_label] = 'public{0}'.format(tree_path.replace('/app/iso/static/public', ''))			    
						    
			    # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
			    # Make histogram Wild Type vs. Variation
			    # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
			    for go_class, go_list in histo_go_pvalue.iteritems():
				label = '{0}_his'.format(go_class)	
				go_l = []
				for g in go_list:
				    wt_p = 0.0
				    var_p = 0.0
				    
				    # Wild Type GO term pvalue
				    if g in wt_go_pvalue[go_class]:
					wt_p = wt_go_pvalue[go_class][g]
					
				    # Variation GO term pvalue
				    if g in var_go_pvalue[go_class]:
					var_p = var_go_pvalue[go_class][g]
					
				    go_l.append({
					'go_term': g,
					'wt': wt_p,
					'var': var_p
				    })
				
				s1 = sorted(go_l, key = itemgetter('var'))
				variation[label] = sorted(s1, key = itemgetter('wt'), reverse = True)
						
			    enrich_json[mirna_id]['variations'][var_type].append(variation)
			    
			    #///////////////////////////////////////////////////////////////
			    # Wild Type and Variation comparison
			    #///////////////////////////////////////////////////////////////
			    comparison = {
				'wt_only_dir_path': 'wt_only', 
				'intersection_dir_path': 'intersection',
				'variation_only_dir_path': 'variation_only'
			    }
			    
			    histo_go_pvalue = { 'bp': [], 'cc': [], 'mf': [] }
			    
			    for c_path, c_label in comparison.iteritems():
				# Get all GO terms for each GO class
				for go_class, go_file in enrich_files.iteritems():
				    f_path = '{0}/{1}'.format(mir[c_path].replace('###', enrich_opt), go_file) 
				    if os.path.isfile(f_path):
					with open(f_path, 'r') as f:
					    d = json.load(f)
					    variation[c_label][go_class] = d    
					    histo_go_pvalue[go_class] = list(set(
						d.keys() + histo_go_pvalue[go_class]))
									    
			    # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
			    # Make histogram Wild Type and Variation comparison
			    # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
			    # For each GO class: bp, cc, mf
			    for go_class, go_list in histo_go_pvalue.iteritems():
				label = 'comp_{0}_his'.format(go_class)		
				
				comparison_go_pvalue = {
				    'wt_only': {},
				    'intersection': {},
				    'variation_only': {},
				}
				
				# For each entity: wt_only, intersection, variation_only
				for c_path, c_label in comparison.iteritems():
				    go_file = go_pvalues_files[go_class]
				    f_path = '{0}/{1}'.format(mir[c_path].replace('###', enrich_opt), go_file)
				    if os.path.isfile(f_path):
					with open(f_path, 'r') as f:
					    d = json.load(f)
					    comparison_go_pvalue[c_label] = d
				go_l = []
				# For each GO term
				for g in go_list:
				    d = {
					'go_term': g,
					'wt_only': 0.0,
					'variation_only': 0.0,
					'intersection': 0.0
				    }
				    
				    # For each entity: wt_only, intersection, variation_only
				    for c_label, go_pvalue in comparison_go_pvalue.iteritems():
					if g in go_pvalue:
					    d[c_label] = go_pvalue[g]
					    
				    go_l.append(d)
				    
				s1 = sorted(go_l, key = itemgetter('wt_only'), reverse = True)
				variation[label] = sorted(s1, key = itemgetter('variation_only'))
	
	    else:
		# Pathway
		for mir in enrich_dirs_list:    
		    mirna_id = mir['wt']['dir_label']
		    wt_path = mir['wt']['enrichment_dir_path'].replace('###', enrich_opt)
		    
		    var_id = mir['variation']['dir_label'].split('__')
		    var_path = mir['variation']['enrichment_dir_path'].replace('###', enrich_opt)
		    
		    wt_pt_pvalue = {}
		    var_pt_pvalue = {}
		    histo_pt_pvalue = []
		    
		    # ==================================================================
		    # ==================================================================
		    # Get all Wild Type and Variation Pathway P-values
		    # ==================================================================
		    # ==================================================================	    
		    f_path = '{0}/enrich_pathway.json'.format(wt_path)
		    if os.path.isfile(f_path):
			with open(f_path, 'r') as f:
			    wt_pt_pvalue = json.load(f)
			    
		    f_path = '{0}/enrich_pathway.json'.format(var_path)
		    if os.path.isfile(f_path):
			with open(f_path, 'r') as f:
			    var_pt_pvalue = json.load(f)
		    
		    
		    # //////////////////////////////////////////////////////////////////
		    # Wild Type
		    # //////////////////////////////////////////////////////////////////
		    if mirna_id not in enrich_json:
			enrich_json[mirna_id] = {
			    'mirna_id': mirna_id,
			    'variations': {
				'substitution': [],
				'isomir': []
			    },
			    'pathways': {},
			    'threshold': math.fabs(math.log10(self.pvalue))
			}
		
		    # Get all Wild Type Pathways
		    f_path = '{0}/enrich.json'.format(wt_path)
		    if os.path.isfile(f_path):
			with open(f_path, 'r') as f:
			    enrich_json[mirna_id]['pathways'] = json.load(f)
			    
		    # //////////////////////////////////////////////////////////////////
		    # Get all variation pathways             
		    # //////////////////////////////////////////////////////////////////
		    if mir['variation']['dir_label'] != '':
			variation = {}
			histo_pt_ids = []
			var_type = None
			matchObj1 = re.match(r'^([0-9]+)_([^_]+)_([^_]+)$', var_id[1], re.M|re.I)
			matchObj2 = re.match(r'^([^_]+)_([^_]+)$', var_id[1], re.M|re.I)
			
			if matchObj1:  
			    var_type = 'substitution'
			    variation['position'] = int(matchObj1.group(1))
			    variation['ref'] = matchObj1.group(2)
			    variation['mod'] = matchObj1.group(3)
			elif matchObj2:
			    var_type = 'isomir'
			    variation['5p'] = matchObj2.group(1)
			    variation['3p'] = matchObj2.group(2)  
		
			if var_type:
			    variation['pathways'] = {}
			    variation['wt_only'] = {}
			    variation['intersection'] = {}
			    variation['variation_only'] = {}  
			    variation['his'] = []
			    variation['comp_his'] = []
			    
			    # Get all Variation GO terms for each GO class
			    f_path = '{0}/enrich.json'.format(var_path) 
			    if os.path.isfile(f_path):
				with open(f_path, 'r') as f:
				    variation['pathways'] = json.load(f) 
				    wt_pt = enrich_json[mirna_id]['pathways'].keys()
				    var_pt = variation['pathways'].keys()
				    histo_pt_ids = list(set(wt_pt + var_pt))
		    
			# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
			# Make histogram Wild Type vs. Variation
			# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
			pt_l = []
			for pathway_id in histo_pt_ids:
			    wt_p = 0.0
			    var_p = 0.0
			    
			    # Wild Type Pathway pvalue
			    if pathway_id in wt_pt_pvalue:
				wt_p = wt_pt_pvalue[pathway_id]
				
			    # Variation Pathway pvalue
			    if pathway_id in var_pt_pvalue:
				var_p = var_pt_pvalue[pathway_id]
				
			    pt_l.append({
				'pathway_id': pathway_id,
				'wt': wt_p,
				'var': var_p
			    })
			    
			s1 = sorted(pt_l, key = itemgetter('var'))
			variation['his'] = sorted(s1, key = itemgetter('wt'), reverse = True)
					    
			enrich_json[mirna_id]['variations'][var_type].append(variation)
			
			#///////////////////////////////////////////////////////////////////
			# Wild Type and Variation comparison
			#///////////////////////////////////////////////////////////////////
			comparison = {
			    'wt_only_dir_path': 'wt_only', 
			    'intersection_dir_path': 'intersection',
			    'variation_only_dir_path': 'variation_only'
			}
			
			histo_pt_pvalue = []
			
			for c_path, c_label in comparison.iteritems():
			    # Get all Pathways
			    f_path = '{0}/enrich.json'.format(mir[c_path].replace('###', enrich_opt))
			    if os.path.isfile(f_path):
				with open(f_path, 'r') as f:
				    variation[c_label] = json.load(f)
				    histo_pt_ids = list(set(variation[c_label].keys() + histo_pt_ids))
				    
			# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
			# Make histogram Wild Type and Variation comparison
			# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
			comparison_pathway_pvalue = {
			    'wt_only': {},
			    'intersection': {},
			    'variation_only': {},
			}
			
			# For each entity: wt_only, intersection, variation_only
			for c_path, c_label in comparison.iteritems():
			    f_path = '{0}/enrich_pathway.json'.format(mir[c_path].replace('###', enrich_opt))
			    if os.path.isfile(f_path):
				with open(f_path, 'r') as f:
				    comparison_pathway_pvalue[c_label] = json.load(f)
				    
			pt_l = []
			for pathway_id in histo_pt_ids:		
			    d = {
				'pathway_id': pathway_id,
				'wt_only': 0.0,
				'variation_only': 0.0,
				'intersection': 0.0
			    }
			    
			    # For each entity: wt_only, intersection, variation_only
			    for c_label, pathway_pvalue in comparison_pathway_pvalue.iteritems():
				if pathway_id in pathway_pvalue:
				    d[c_label] = pathway_pvalue[pathway_id]
				    
			    pt_l.append(d)
			    
			s1 = sorted(pt_l, key = itemgetter('wt_only'), reverse = True)
			variation['comp_his'] = sorted(s1, key = itemgetter('variation_only'))		
		
	    master_enrich_json[enrich_opt] = enrich_json
	    
	with open(en_path, 'w') as f:
	    f.write(json.dumps(master_enrich_json, indent=4))  
	    
	    
    def ePrint(self, message):
        print >> sys.stderr, message
        
        
    def ePrintClose(self, message):
        print >> sys.stderr, message  
        enrichment_jsonsys.exit(1)
        
        
     