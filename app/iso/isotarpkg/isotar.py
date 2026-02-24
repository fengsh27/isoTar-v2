#!/usr/bin/env python 
# -*- coding: utf-8 -*-
"""
================================================================================
    Author: Rosario Distefano
    Email: rosario.distefano.ict@gmail.com
    Project: isoTar - isoform Targeting
================================================================================
"""

''' Libraries '''
import os.path
import cmdargs as cmda
import pred_analysis as pr
import enrich_analysis as en
import sys

class Isotar():
    def __init__(self, mirna_path, target_path, output_dir, map_path, go_path,
                 pre_mature_path, pre_mat_acc_path, nm_to_geneid_path,
                 obo_dag_path, gene_2_go_path, pathways_path, gene_2_pathway_path):
	
        ''' isoTar constructor '''
        # Resource paths
        self.mirna_path = mirna_path
        self.target_path = target_path
        self.output_dir = output_dir
        self.go_path = go_path
        self.map_path = map_path 
        self.pre_mature_path = pre_mature_path
	self.pre_mat_acc_path = pre_mat_acc_path
	self.nm_to_geneid_path = nm_to_geneid_path
	self.obo_dag_path = obo_dag_path
	self.gene_2_go_path = gene_2_go_path
	self.pathways_path = pathways_path
	self.gene_2_pathway_path = gene_2_pathway_path
        
        if not os.path.isfile(self.mirna_path):
            self.ePrintClose("Error. No mirnas.fa file found.")
        
        if not os.path.isfile(self.target_path):
            self.ePrintClose("Error. No 3utr.fa file found.")
        
        if not os.path.isdir(self.output_dir):
            self.ePrintClose("Error. No output directory found.")
        
        if not os.path.isfile(self.map_path):
            self.ePrintClose("Error. No targets_go_map.tab file found.")
            
        if not os.path.isfile(self.go_path):
            self.ePrintClose("Error. No go.json file found.")
            
        if not os.path.isfile(self.pre_mature_path):
            self.ePrintClose("Error. No pre_mature_path.json file found.")
	    
	if not os.path.isfile(self.pre_mat_acc_path):
	    self.ePrintClose("Error. No pre_mat_acc.json file found.")
	
	if not os.path.isfile(self.nm_to_geneid_path):
	    self.ePrintClose("Error. No nm_to_geneid.tab file found.")
	    
	if not os.path.isfile(self.obo_dag_path):
	    self.ePrintClose("Error. No go-basic.obo file found.")
	
	if not os.path.isfile(self.gene_2_go_path):
	    self.ePrintClose("Error. No gene2go file found.")
	
	if not os.path.isfile(self.pathways_path):
	    self.ePrintClose("Error. No pathway.json file found.")
	
	if not os.path.isfile(self.gene_2_pathway_path):
	    self.ePrintClose("Error. No genename_pathways.json file found.")
	    	    
        # Define the Arguments Parser
        self.cmd_l = cmda.CommandLineArgsParser()        
    
    
    def ePrint(self, message):
        print >> sys.stderr, message
    
        
    def ePrintClose(self, message):
        print >> sys.stderr, message  
        sys.exit(1)
        
        
    def run(self):
        ''' isoTar run '''
        # Prediction analysis
        if self.cmd_l.prediction == 'yes':
            pred_analysis = pr.PredictionAnalysis(self.mirna_path,
                                                  self.target_path,
                                                  self.output_dir,
                                                  self.pre_mature_path,
                                                  self.cmd_l,
	                                          self.pre_mat_acc_path,
	                                          self.nm_to_geneid_path)
        # Enrichment analysis
        if self.cmd_l.enrichment == 'yes':
            enrich_analysis = en.EnrichmentAnalysis(self.output_dir,
                                                    self.nm_to_geneid_path,
                                                    self.cmd_l,
	                                            self.obo_dag_path,
	                                            self.gene_2_go_path,
	                                            self.go_path,
	                                            self.pathways_path,
	                                            self.gene_2_pathway_path)
      	
        