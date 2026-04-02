# -*- coding: utf-8 -*-

''' Libraries '''
import mimerender
import re
import datetime
from pymongo import MongoClient
import pymongo
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
import multiprocessing as mp
import math
from operator import itemgetter
import csv

class TargetsProcessing():
        
    def __init__(self, path, cores):
        if not os.path.isfile(path) or not os.path.isfile(path): 
            print 'Error. Please specify the Target sequence(s) [ FASTA format ] file.'
            sys.exit(1)
        
        if not isinstance(cores, int) or int(cores) <= 0:
            print 'Error. Please use an integer to specify the number of cores'
            sys.exit(1)      
            
        self.debug = False 
          
        # Get the file basename and extension
        self.basename, self.extension = os.path.splitext(os.path.basename(path))        
          
        self.path = path
        self.targets_list = []
        self.targets_split_list = []
        self.no_targets = 0
        self.no_targets_mirza_g = 0
        self.longest_target_seq = 0
        self.no_cores = 0
        self.target_split_details = []
        self.target_mirza_g_split_details = []
        self.targets_targetscan_split_details = []
        self.targets_map = []
        self.targets_enst_to_refseq_map = {}
        self.no_cores = cores
        self.targets_map_gene_id_to_refseq = {}
        
        
    def setDebugging(self, debug):
        self.debug = debug    
        
        
    def deleteTargetsSupportFile(self):
        
        # Remove all target splitted files
        for f in self.target_split_details:
            if os.path.isfile(f['path']):
                os.remove(f['path'])  # remove the file     
                
        # Remove all target splitted files
        for f in self.target_mirza_g_split_details:
            if os.path.isfile(f['path']):
                os.remove(f['path'])  # remove the file 
            
        
    def targetsEnstToRefseqMap(self, path):    
        targets_map = {}
        if os.path.isfile(path):
            with open(path, 'r') as f:
                handler = csv.reader(f, delimiter='\t')
                next(handler, None)
                for line in handler:
                    if line[0].strip() != '' and line[1].strip() != '':
                        targets_map[line[0].strip()] = line[1].strip()

        self.targets_enst_to_refseq_map = targets_map        
    
    
    def prepareTargetscanTargets(self, utr_path, bls_bin_path):    
        
        tmp_targets = []
        # No splitted target
        if os.path.isdir(utr_path):
            tmp_targets = [ {
                'path': '%s/%s' % (utr_path, f), 
                'bls_path': '%s/targetscan_median_bls_bins_part_%s.txt' % (bls_bin_path,
                                                                           re.match(r'targetscan_utr_part_([0-9]+)\.txt', f).group(1)),
                'idx': int(re.match(r'targetscan_utr_part_([0-9]+)\.txt', f).group(1))
            } for f in os.listdir(utr_path) if re.match(r'targetscan_utr_part_[0-9]+\.txt', f)]

        targets = sorted(tmp_targets, 
                         key=itemgetter('idx'), 
                         reverse=False)   

        self.targets_targetscan_split_details = targets
        
    
    def targetsMartMap(self, path):    
        targets_map = {}
        if os.path.isfile(path):
            with open(path, 'r') as f:
                handler = csv.reader(f, delimiter='\t')
                # Avoid the header
                next(handler, None)
                for line in handler:
                    row = {
                        'refseq_id': line[0].strip(),                        
                        'gene_id': line[1].strip(),
                        'gene_name': line[2].strip()
                    }
                    
                    if row['refseq_id'] != '':
                        targets_map[row['refseq_id']] = row
                            
        return targets_map
    
    
    def searchIntoTargetsMap(self, field_name, value_to_find, targets_list):
        
        row = next((item for item in targets_list if value_to_find in item[field_name]), {})
        return row
    
        
    def getTargets(self):
        
        self.targets_map = self.targetsMartMap('/app/iso/static/public/resources/nm_to_geneid.tab')
        
        content = []
        tmp_targets_list = []
        
        # Keep all identified identifier ID
        tmp_targets_identifier_list = []
        
        # Read the targets file
        with open(self.path, 'r') as f:
            content = f.readlines()
            
        for i in xrange(0, len(content)):
            # Get the target header
            if re.match(r'^>.*$', content[i], re.M|re.I):
                matchObj = re.match(r'^>.*([NM]{2}[^\s]+)\s.*$', content[i], re.M)
                row = {}
                gene_id = ''
                gene_name = ''
                if matchObj:
                    refseq = matchObj.group(1)
                    row = self.targets_map[refseq] if refseq in self.targets_map else {}
                    if row :
                        # Gene ID
                        if 'gene_id' in row:
                            gene_id = row['gene_id']
                            
                            if gene_id not in self.targets_map_gene_id_to_refseq:
                                self.targets_map_gene_id_to_refseq[gene_id] = []
                            
                            if refseq not in self.targets_map_gene_id_to_refseq[gene_id]:
                                self.targets_map_gene_id_to_refseq[gene_id].append(refseq)                                                                
                        # Gene name
                        if 'gene_name' in row:
                            gene_name = row['gene_name']
                                            
                    target = { 
                        'header': content[i].replace('\n',''), 
                        'identifier': refseq,
                        'seq': '',
                        'length': 0,
                        'gene_id': gene_id,
                        'gene_name': gene_name
                    }
                    i += 1
                    tmp_seq = ''
                    
                    # Get the target sequence
                    while i < len(content) and \
                          not re.match(r'^>.*$', content[i], re.M|re.I):
                        tmp_seq += content[i].replace('\n','')
                        i += 1
                        
                    i -= 1
                    target['seq'] = tmp_seq
                    target['length'] = len(tmp_seq)
                    
                    # Get the longest target sequence
                    if self.longest_target_seq < len(tmp_seq):
                        self.longest_target_seq = len(tmp_seq)
                    
                    # Store the target details
                    if target['identifier'] != '' and \
                       target['identifier'] not in tmp_targets_identifier_list:
                        tmp_targets_identifier_list.append(target['identifier'])
                        tmp_targets_list.append(target)
        
        # Decending sorting by sequence length
        self.targets_list = sorted(tmp_targets_list, 
                                   key=itemgetter('length'), 
                                   reverse=True)
    
        # No. of targets
        self.no_targets = len(self.targets_list)
                    
    def printTargetsDetails(self):
        for t in self.targets_list:
            print t
            
            
    def prepareTargets(self):
        
        if self.no_cores > 1:
            # Main tools
            tmp_targets_split_list = [[] for i in xrange(0, self.no_cores)]
            
            for i in xrange(0, self.no_cores):
                for j in xrange(0, len(self.targets_list)):
                    idx = (j*self.no_cores) + i
                    if idx < len(self.targets_list):
                        tmp_targets_split_list[i].append(self.targets_list[idx])
            
            self.targets_split_list = tmp_targets_split_list      
            
            #basic_path = self.path[:-len(self.extension)]
            basic_path = '/input/3utr'
            
            # Statistics
            for i in xrange(0, len(self.targets_split_list)):
                
                if self.debug:
                    print 'Chuck %d, size: %d' % (i, 
                                                  len(self.targets_split_list[i]))
                
                if len(self.targets_split_list[i]) > 0:
                    
                    path = '%s_part_%d%s' % (basic_path, i, self.extension)
                    
                    self.target_split_details.append({
                        'path': path,
                        'longest_target_seq': len(self.targets_split_list[i][0]['seq'])
                    })
                    
                    with open(path, 'w') as t:            
                        for j in xrange(0, len(self.targets_split_list[i])):
                            line = '%s\n' % self.targets_split_list[i][j]['header']
                            t.write(line)
                            line = '%s\n' % self.targets_split_list[i][j]['seq']
                            t.write(line)    
            
            if self.debug:
                print self.target_split_details