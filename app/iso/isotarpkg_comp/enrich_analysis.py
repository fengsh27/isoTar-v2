# -*- coding: utf-8 -*-

''' Libraries '''
import re
import argparse
import os, errno
import os.path
import sys
import psutil
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
        self.output_enrichment_dir = '{0}/enrichment_comp'.format(output_dir)
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
        
        # Delete the enrichment directory
        utl.removeDirectory(self.output_enrichment_dir)
	utl.createDirectory(self.output_enrichment_dir)
	
	# Create all specified enrichment type directories
	for enrich_opt in self.enrich_type:
	    utl.createDirectory('{0}/{1}'.format(
	    self.output_enrichment_dir, enrich_opt))	
	
	enrich_dirs_list = {}
	
	for enrich_opt in self.enrich_type:
		
	    # ==============================================================
	    # GO term Enrichment Analysis
	    # ==============================================================
	    if 'go' == enrich_opt:
		obo_dag = GODag(self.obo_dag_path)
		geneid2gos_human = read_ncbi_gene2go(self.gene_2_go_path, taxids=[9606])
		
		# Set the GO term object
		GO = go.GO(self.nm_to_geneid_path, 
			   cmd_l, 
			   self.go_path,
			   '{0}/{1}'.format(self.output_enrichment_dir, enrich_opt))	
		    
		# Workers queue
		go_queue = mp.Queue()
			    
		# Execute the enrichment analysis for all miRNA(s) 
		GO.executeEnrichmentAnalysis(cmd_l.no_cpu_cores, 
			                     go_queue,
			                     obo_dag, 
			                     geneid2gos_human,
			                     utl)
		   
		#self.prepareEnrichedJsonFile(GO.support)
		enrich_dirs_list[enrich_opt] = GO.support
	    else:
		# ==========================================================
		# Pathway Enrichment Analysis
		# ==========================================================		    
		# Set the Pathway object
		Pathway = pt.Pathway(self.nm_to_geneid_path, enrich_opt, 
	                             cmd_l, self.pathways_path,
	                             self.gene_2_pathway_path, 
		                     '{0}/{1}'.format(self.output_enrichment_dir, enrich_opt))
		
		# Workers queue
		pt_queue = mp.Queue()
		
		# Execute the pathway enrichment analysis for all miRNA(s) 
		Pathway.executeEnrichmentAnalysis(cmd_l.no_cpu_cores, 
	                                          pt_queue,
	                                          Pathway.worker,
	                                          enrich_opt,
		                                  utl)
		
		enrich_dirs_list[enrich_opt] = Pathway.support
		
	self.prepareEnrichedJsonFile(enrich_dirs_list)		
                
                
    def ePrint(self, message):
        print >> sys.stderr, message
        
        
    def ePrintClose(self, message):
        print >> sys.stderr, message  
        enrichment_jsonsys.exit(1)
        
        
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
	
		list1_go_pvalue = { 'bp': [], 'cc': [], 'mf': [] }
		list2_go_pvalue = { 'bp': [], 'cc': [], 'mf': [] }	
		histo_go_pvalue = { 'bp': [], 'cc': [], 'mf': [] }
		
		'''
		////////////////////////////////////////////////////////////////////////
		Check list1
		////////////////////////////////////////////////////////////////////////
		'''
		if 'enrichment_list1_dir_path' in enrich_dirs_list[enrich_opt]:
		    dir_path = enrich_dirs_list[enrich_opt]['enrichment_list1_dir_path']
			
		    # Get all list1 GO term pvalues
		    for go_class, go_file in go_pvalues_files.iteritems():
			f_path = '{0}/{1}'.format(dir_path, go_file)
			if os.path.isfile(f_path):
			    with open(f_path, 'r') as f:
				d = json.load(f)
				list1_go_pvalue[go_class] = d	
				
		    enrich_json['list1'] = {
			'entity': enrich_dirs_list[enrich_opt]['name_list1'],
			'bp': {},
			'cc': {},
			'mf': {},
			'threshold': math.fabs(math.log10(self.pvalue)),
			'bp_tree': '',
			'cc_tree': '',
			'mf_tree': ''
		    }
		    
		    # Get all list1 GO terms for each GO class
		    for go_class, go_file in enrich_files.iteritems():
			f_path = '{0}/{1}'.format(dir_path, go_file) 
			if os.path.isfile(f_path):
			    with open(f_path, 'r') as f:
				d = json.load(f)
				enrich_json['list1'][go_class] = d
				# List1 GO Tree
				tree_path = '{0}/{1}_go_tree_{2}.pdf'.format(
				    dir_path, 
				    enrich_dirs_list[enrich_opt]['name_list1'],
				    go_class) 
				tree_label = '{0}_tree'.format(go_class)
				if os.path.isfile(tree_path):
				    enrich_json['list1'][tree_label] = 'public{0}'.format(
				        tree_path.replace('/app/iso/static/public', ''))
				    
		
		'''
		////////////////////////////////////////////////////////////////////////
		Check list2
		////////////////////////////////////////////////////////////////////////
		'''
		if 'enrichment_list2_dir_path' in enrich_dirs_list[enrich_opt]:
		    dir_path = enrich_dirs_list[enrich_opt]['enrichment_list2_dir_path']
			
		    # Get all list2 GO term pvalues
		    for go_class, go_file in go_pvalues_files.iteritems():
			f_path = '{0}/{1}'.format(dir_path, go_file)
			if os.path.isfile(f_path):
			    with open(f_path, 'r') as f:
				d = json.load(f)
				list2_go_pvalue[go_class] = d	
				
		    enrich_json['list2'] = {
			'entity': enrich_dirs_list[enrich_opt]['name_list2'],
			'bp': {},
			'cc': {},
			'mf': {},
			'threshold': math.fabs(math.log10(self.pvalue)),
			'bp_tree': '',
			'cc_tree': '',
			'mf_tree': '',
			'wt_only': {},
			'intersection': {},
			'variation_only': {},
			'bp_his': [],
			'cc_his': [],
			'mf_his': [],
			'comp_bp_his': [],
			'comp_cc_his': [],
			'comp_mf_his': [],
			'bp_wt_vs_variation_tree': '',
			'cc_wt_vs_variation_tree': '',
			'mf_wt_vs_variation_tree': ''	        
		    }
		    
		    # Get all list2 GO terms for each GO class
		    for go_class, go_file in enrich_files.iteritems():
			f_path = '{0}/{1}'.format(dir_path, go_file) 
			if os.path.isfile(f_path):
			    with open(f_path, 'r') as f:
				d = json.load(f)
				enrich_json['list2'][go_class] = d
				list1_go = enrich_json['list1'][go_class].keys()
				histo_go_pvalue[go_class] = list(set(list1_go + d.keys()))
				
				# List2 GO Tree
				tree_path = '{0}/{1}_go_tree_{2}.pdf'.format(
				    dir_path, 
				    enrich_dirs_list[enrich_opt]['name_list2'],
				    go_class) 
				tree_label = '{0}_tree'.format(go_class)
				if os.path.isfile(tree_path):
				    enrich_json['list2'][tree_label] = 'public{0}'.format(tree_path.replace('/app/iso/static/public', ''))   
				
				# Comparison GO Tree
				tree_path = '{0}/go_tree_comp_{1}.pdf'.format(
				    enrich_dirs_list[enrich_opt]['enrichment_comp_dir_path'], 
				    go_class) 
				tree_label = '{0}_wt_vs_variation_tree'.format(go_class)
				if os.path.isfile(tree_path):
				    enrich_json['list2'][tree_label] = 'public{0}'.format(tree_path.replace('/app/iso/static/public', ''))
	    
		    # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
		    # Make histogram list1 vs. list2
		    # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
		    for go_class, go_list in histo_go_pvalue.iteritems():
			label = '{0}_his'.format(go_class)	
			go_l = []
			for g in go_list:
			    list1_p = 0.0
			    list2_p = 0.0
			    
			    # List1 GO term pvalue
			    if g in list1_go_pvalue[go_class]:
				list1_p = list1_go_pvalue[go_class][g]
				
			    # List2 GO term pvalue
			    if g in list2_go_pvalue[go_class]:
				list2_p = list2_go_pvalue[go_class][g]
				
			    go_l.append({
				'go_term': g,
				'wt': list1_p,
				'var': list2_p
			    })
			
			s1 = sorted(go_l, key = itemgetter('var'))
			enrich_json['list2'][label] = sorted(s1, key = itemgetter('wt'), reverse = True)
					    
		    #///////////////////////////////////////////////////////////////////
		    # List1 and List2 comparison
		    #///////////////////////////////////////////////////////////////////
		    
		    if 'enrichment_comp_dir_path' in enrich_dirs_list[enrich_opt]:
			comparison = {
			    'enrichment_comp_list1_only_dir_path': 'wt_only', 
			    'enrichment_comp_intersection_dir_path': 'intersection',
			    'enrichment_comp_list2_only_dir_path': 'variation_only'
			}
			
			histo_go_pvalue = { 'bp': [], 'cc': [], 'mf': [] }
			
			for c_path, c_label in comparison.iteritems():
			    # Get all GO terms for each GO class
			    for go_class, go_file in enrich_files.iteritems():
				f_path = '{0}/{1}'.format(enrich_dirs_list[enrich_opt][c_path], go_file) 
				if os.path.isfile(f_path):
				    with open(f_path, 'r') as f:
					d = json.load(f)
					enrich_json['list2'][c_label][go_class] = d    
					histo_go_pvalue[go_class] = list(set(
					    d.keys() + histo_go_pvalue[go_class]))
									
			# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
			# Make histogram list1 and list2 comparison
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
				f_path = '{0}/{1}'.format(enrich_dirs_list[enrich_opt][c_path], go_file)
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
			    enrich_json['list2'][label] = sorted(s1, key = itemgetter('variation_only'))
			    
	    else:
		# Pathway
		list1_data_pvalue = {}
		list2_data_pvalue = {}	
		histo_data_pvalue = []
		
		'''
		////////////////////////////////////////////////////////////////////////
		Check list1
		////////////////////////////////////////////////////////////////////////
		'''
		if 'enrichment_list1_dir_path' in enrich_dirs_list[enrich_opt]:
		    dir_path = enrich_dirs_list[enrich_opt]['enrichment_list1_dir_path']
		    
		    # Get all list1 Pathway pvalues
		    f_path = '{0}/enrich_pathway.json'.format(dir_path)
		    if os.path.isfile(f_path):
			with open(f_path, 'r') as f:
			    d = json.load(f)
			    list1_data_pvalue = d	
				
		    enrich_json['list1'] = {
			'entity': enrich_dirs_list[enrich_opt]['name_list1'],
			'pathways': {},
			'threshold': math.fabs(math.log10(self.pvalue)),
		    }
		    
		    # Get all list1 Pathway details
		    f_path = '{0}/enrich.json'.format(dir_path) 
		    if os.path.isfile(f_path):
			with open(f_path, 'r') as f:
			    d = json.load(f)
			    enrich_json['list1']['pathways'] = d
			    
		'''
		////////////////////////////////////////////////////////////////////////
		Check list2
		////////////////////////////////////////////////////////////////////////
		'''
		if 'enrichment_list2_dir_path' in enrich_dirs_list[enrich_opt]:
		    dir_path = enrich_dirs_list[enrich_opt]['enrichment_list2_dir_path']
			
		    # Get all list2 Pathway pvalues
		    f_path = '{0}/enrich_pathway.json'.format(dir_path)
		    if os.path.isfile(f_path):
			with open(f_path, 'r') as f:
			    d = json.load(f)
			    list2_data_pvalue = d	
				
		    enrich_json['list2'] = {
			'entity': enrich_dirs_list[enrich_opt]['name_list2'],
			'pathways': {},
			'threshold': math.fabs(math.log10(self.pvalue)),
			'wt_only': {},
			'intersection': {},
			'variation_only': {},
			'his': [],
			'comp_his': []
		    }	    
		    
		    # Get all list2 Pathway details
		    f_path = '{0}/enrich.json'.format(dir_path) 
		    if os.path.isfile(f_path):
			with open(f_path, 'r') as f:
			    enrich_json['list2']['pathways'] = json.load(f)
			    list1_pathway_ids = enrich_json['list1']['pathways'].keys()
			    list2_pathway_ids = enrich_json['list2']['pathways'].keys()
			    histo_data_pvalue = list(set(list1_pathway_ids + list2_pathway_ids))			    
		    
		    
		
		    # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
		    # Make histogram list1 vs. list2
		    # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
		    pt_l = []
		    
		    # For each pathway
		    for pathway_id in histo_data_pvalue:
			list1_p = 0.0
			list2_p = 0.0
			
			# List1 Pathway pvalue
			if pathway_id in list1_data_pvalue:
			    list1_p = list1_data_pvalue[pathway_id]
			    
			# List2 Pathway pvalue
			if pathway_id in list2_data_pvalue:
			    list2_p = list2_data_pvalue[pathway_id]
			    
			pt_l.append({
			    'pathway_id': pathway_id,
			    'wt': list1_p,
			    'var': list2_p
			})
		    
		    s1 = sorted(pt_l, key = itemgetter('var'))
		    enrich_json['list2']['his'] = sorted(s1, key = itemgetter('wt'), reverse = True)
		    
		    
		    #///////////////////////////////////////////////////////////////////
		    # List1 and List2 comparison
		    #///////////////////////////////////////////////////////////////////
		    
		    if 'enrichment_comp_dir_path' in enrich_dirs_list[enrich_opt]:
			comparison = {
			    'enrichment_comp_list1_only_dir_path': 'wt_only', 
			    'enrichment_comp_intersection_dir_path': 'intersection',
			    'enrichment_comp_list2_only_dir_path': 'variation_only'
			}
			
			histo_data_pvalue = []
			
			for c_path, c_label in comparison.iteritems():
			    # Get all Pathway details for ech subclass (e.g.: wt_only, ..)
			    f_path = '{0}/enrich.json'.format(enrich_dirs_list[enrich_opt][c_path]) 
			    if os.path.isfile(f_path):
				with open(f_path, 'r') as f:
				    enrich_json['list2'][c_label] = json.load(f)
				    pathway_ids = enrich_json['list2'][c_label].keys()    
				    histo_data_pvalue = list(set(pathway_ids + histo_data_pvalue))
									
			# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
			# Make histogram list1 and list2 comparison
			# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
			comparison_data_pvalue = {
			    'wt_only': {},
			    'intersection': {},
			    'variation_only': {},
			}
			
			# For each entity: wt_only, intersection, variation_only
			for c_path, c_label in comparison.iteritems():
			    f_path = '{0}/enrich_pathway.json'.format(enrich_dirs_list[enrich_opt][c_path])
			    if os.path.isfile(f_path):
				with open(f_path, 'r') as f:
				    d = json.load(f)
				    comparison_data_pvalue[c_label] = d
			pt_l = []
			# For each pathway
			for pathway_id in histo_data_pvalue:
			    d = {
				'pathway_id': pathway_id,
				'wt_only': 0.0,
				'variation_only': 0.0,
				'intersection': 0.0
			    }
			    
			    # For each entity: wt_only, intersection, variation_only
			    for c_label, pathway_pvalue in comparison_data_pvalue.iteritems():
				if pathway_id in pathway_pvalue:
				    d[c_label] = pathway_pvalue[pathway_id]
				    
			    pt_l.append(d)
			    
			s1 = sorted(pt_l, key = itemgetter('wt_only'), reverse = True)
			enrich_json['list2']['comp_his'] = sorted(s1, key = itemgetter('variation_only'))
			    
	    master_enrich_json[enrich_opt] = enrich_json
	    
	with open(en_path, 'w') as f:
	    f.write(json.dumps(master_enrich_json, indent=4))