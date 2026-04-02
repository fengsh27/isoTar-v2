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
import mirmap
import csv
from operator import itemgetter
import gzip

class ITUtils():
        
    def __init__(self, cores):
        self.debug = False
        self.no_cores = cores
    
    def ePrint(self, message):
        print >> sys.stderr, message
        
    def ePrintClose(self, message):
        print >> sys.stderr, message  
        sys.exit(1)
	
    def setDebugging(self, debug):
        self.debug = debug
	            
    def removeDirectory(self, path):
        """ param <path> could either be relative or absolute. """
        if os.path.exists(path):
            if os.path.isfile(path):
                os.remove(path)  # remove the file
            elif os.path.isdir(path):
                try:
                    shutil.rmtree(path) # remove dir and all contains
                except OSError as e:
                    pass         
            else:
                pass

    def removeDirectoryContent(self, path):
        """ param <path> could either be relative or absolute. """
        if os.path.exists(path):
	    for p in os.listdir(path):
		cpath = os.path.join(path, p)
                if os.path.isfile(cpath):
                    os.remove(cpath)  # remove the file
                elif os.path.isdir(cpath):
                    try:
                        shutil.rmtree(cpath) # remove dir and all contains
                    except OSError as e:
                        if e.errno != errno.EEXIST:
                            raise         
                else:
                    raise ValueError("file {} is not a file or dir.".format(cpath))
        
    def createDirectory(self, path):
        # Create the directory
        try:
            os.makedirs(path)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise         
    
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
    
    def mirmapRuns(self, mirna_seq, target_seq):
        
        # Check U-T nucleotides
        target_u = target_seq.find('U')
        mirna_u = mirna_seq.find('U')
                
        # Replace miRNA's U with T
        if target_u == -1 and mirna_u != -1:
            mirna_seq = mirna_seq.replace('U', 'T')
            
        # Replace miRNA's T with U
        if target_u != -1 and mirna_u == -1:
            mirna_seq = mirna_seq.replace('T', 'U')
            
        output = ''
        
        try:
            mim = mirmap.mm(target_seq, mirna_seq)
            mim.find_potential_targets_with_seed(
                allowed_lengths=[7,8], allowed_gu_wobbles={7:0,8:0},
                allowed_mismatches={7:0,8:0}, take_best=True)
            mim.end_sites                                    # Coordinate(s) (3' end) of the target site on the target sequence
            mim.eval_tgs_au(with_correction=False)           # TargetScan features manually evaluated with
            mim.eval_tgs_pairing3p(with_correction=False)    # a non-default parameter.
            mim.eval_tgs_position(with_correction=False)
            mim.eval_prob_binomial()                                # mim's attribute: the feature is automatically computed
            output = mim.report()  
        except ValueError as err:
            pass 
        
        return output
        
    
    
    def worker(self, command, idx, tool, mirna_obj, queue):
        """thread worker function"""
        if self.debug:
            print 'Worker No. %d. is executing %s on miRNA %s.' % (idx, 
                                                                   tool['tool'],
                                                                   mirna_obj['identifier'])
        pre_time = datetime.datetime.now()
        command_status = self.executeCommandLine(command)
        next_time = datetime.datetime.now()
        
        delta_time = next_time - pre_time
        running_time = int(delta_time.total_seconds() * 10**6)
        
        if self.debug:
            if command_status != -1:
                print '\t*** %s has been executed on miRNA %s (Worker No. %d, Time \
                (microseconds): %d)' % (tool['tool'], 
                                        mirna_obj['identifier'], 
                                        idx,
                                        running_time)
            else:
                print '\t*** An error occurred during %s execution \
                on miRNA %s (Worker No. %d)' %  (tool['tool'], 
                                                 mirna_obj['identifier'], 
                                                 idx) 
                
        queue.put(({ 
            'tool': tool['tool'], 
            'miRNA': mirna_obj['identifier'], 
            'time': running_time,
            'error': ''
        }, idx))
            
        return
    
    
    def mergeWorker(self, idx, dirs_list):
        """thread worker function"""
        
        # For each output directory
        for d in dirs_list:
	    # Final output file 
	    o_file = 'output.txt'
	    suffix = '.txt'
	    n_header = 1
	    n_tail = 0
	    if d['tool'] == 'PITA':
		o_file = 'output_pita_results.tab'
		suffix = '_pita_results.tab'
		
	    if d['tool'] in ['RNAhybrid', 'miRmap']:
		n_header = 0
		
	    if d['tool'] == 'miRanda':
		n_header = 32
		n_tail = 2
		
	    m_path = '{0}/{1}'.format(d['dir_path'], o_file)
	    
	    # Prepare the final output file merging all output parts
	    with open(m_path, 'w') as m_out:
		header = False
		for f_idx in xrange(0, d['n_files']):
		    f_path = '{0}/output_part_{1}{2}'.format(d['dir_path'], f_idx, suffix)
		    
		    if os.path.isfile(f_path):			
			with open(f_path, 'r') as f:
			    lines = f.readlines()
			    if len(lines) > 0:
				if not header:
				    header = True
				    m_out.write(''.join(lines[: len(lines)-n_tail]))
				else:
				    m_out.write(''.join(lines[n_header: len(lines)-n_tail]))
				    
			# Delete the i-th output file
			self.removeDirectory(f_path)
			if d['tool'] == 'PITA':
			    d_path = '{0}/output_part_{1}'.format(d['dir_path'], f_idx)
			    self.removeDirectory(d_path)
				    
	    
	    if d['tool'] == 'PITA':
		o_file = 'output_pita_results_targets.tab'
		m_path = '{0}/{1}'.format(d['dir_path'], o_file)
		suffix = '_pita_results_targets.tab'
		
		# Prepare the final output file merging all output parts
		with open(m_path, 'w') as m_out:
		    header = False
		    for f_idx in xrange(0, d['n_files']):
			f_path = '{0}/output_part_{1}{2}'.format(d['dir_path'], f_idx, suffix)
			
			if os.path.isfile(f_path):			
			    with open(f_path, 'r') as f:
				lines = f.readlines()
				if len(lines) > 0:
				    if not header:
					header = True
					m_out.write(''.join(lines[: len(lines)-n_tail]))
				    else:
					m_out.write(''.join(lines[n_header: len(lines)-n_tail]))
			    # Delete the i-th output file
			    self.removeDirectory(f_path)
			    
	    if d['tool'] == 'TargetScan':
		o_file = 'output_sort.txt'
		m_path = '{0}/{1}'.format(d['dir_path'], o_file)
		
		# Prepare the final output file merging all output parts
		with open(m_path, 'w') as m_out:
		    header = False
		    for f_idx in xrange(0, d['n_files']):
			f_path = '{0}/output_sort_part_{1}{2}'.format(d['dir_path'], f_idx, suffix)
			
			if os.path.isfile(f_path):			
			    with open(f_path, 'r') as f:
				lines = f.readlines()
				if len(lines) > 0:
				    if not header:
					header = True
					m_out.write(''.join(lines[: len(lines)-n_tail]))
				    else:
					m_out.write(''.join(lines[n_header: len(lines)-n_tail]))
			    # Delete the i-th output file
			    self.removeDirectory(f_path)
		
        return    
    
    
    def TargetScanWorker(self, tool, idx, mirna_obj, targets_list, queue, exec_type):
        """thread worker function"""
        if self.debug:
            print 'Worker No. %d. is executing %s on miRNA %s.' % (idx, 
                                                                   tool['tool'],
                                                                   mirna_obj['identifier'])
        
        pre_time_glob = datetime.datetime.now()
        for target in targets_list:            
            pre_time = datetime.datetime.now()
            # Target file path
            target_path = target['path']   
            target_bls_path = target['bls_path']
            
            # Output file path
            output_f_path = ''
            output2_f_path = ''
            
            if exec_type == 'single':
                output_f_path = tool['output_filename']
                output2_f_path = tool['output_filename2']                
            else:
                output_f_path = tool['output_splitted_filename'] % target['idx'] 
                output2_f_path = tool['output_splitted_filename2'] % target['idx']
            
            # miRNA FASTA file path
            miRNA_f_path = '%s/mirna_targetscan.txt' % mirna_obj['output_path']
            
            miRNA_tool_dir_path = '%s/%s' % (mirna_obj['output_path'], 
	                                     tool['tool'])    
            
            targetscan_output = '%s/%s' % (miRNA_tool_dir_path,
                                           output_f_path)
            targetscan_output2 = '%s/%s' % (miRNA_tool_dir_path,
                                           output2_f_path)            
            # Execute targetscan_70.pl
            command = tool['cmd'] % (miRNA_f_path, 
                                     target_path,
                                     targetscan_output)  
            command_status = self.executeCommandLine(command)
            
            next_time = datetime.datetime.now()
            delta_time = next_time - pre_time
            running_time = int(delta_time.total_seconds() * 10**6)
        
            if self.debug:
                if command_status != -1:
                    print '\t*** %s has been executed on miRNA %s (Worker No. %d, Time \
                    (microseconds): %d)' % (tool['tool'], 
                                            mirna_obj['identifier'], 
                                            idx,
                                            running_time)
                else:
                    print '\t*** An error occurred during %s execution (CMD 1) \
                    on miRNA %s (Worker No. %d)' %  (tool['tool'], 
                                                     mirna_obj['identifier'], 
                                                     idx)    
            # targetscan_70_BL_PCT.pl      
            if os.path.isfile(targetscan_output):
                pre_time = datetime.datetime.now()
                command = tool['cmd2'] % (miRNA_f_path,
                                          targetscan_output, 
                                          target_bls_path,
                                          targetscan_output2) 
                
                command_status = self.executeCommandLine(command)
                next_time = datetime.datetime.now()
                delta_time = next_time - pre_time
                running_time += int(delta_time.total_seconds() * 10**6) 
                
                if self.debug:
                    if command_status != -1:
                        print '\t*** %s has been executed on miRNA %s (Worker No. %d, Time \
                        (microseconds): %d)' % (tool['tool'], 
                                                mirna_obj['identifier'], 
                                                idx,
                                                running_time)
                    else:
                        print '\t*** An error occurred during %s execution (CMD 2)\
                        on miRNA %s (Worker No. %d)' %  (tool['tool'], 
                                                         mirna_obj['identifier'], 
                                                         idx)                
        
        next_time_glob = datetime.datetime.now()
        delta_time_glob = next_time_glob - pre_time_glob
        running_time_glob = int(delta_time_glob.total_seconds() * 10**6)        
            
        queue.put(({
            'tool': tool['tool'], 
            'miRNA': mirna_obj['identifier'], 
            'time': running_time_glob,
            'error': ''
        }, idx))
            
        return     
    
    
    def miRmapWorker(self, idx, tool, mir_obj, targets_list, output_path, queue):
        """thread worker function"""
        if self.debug:
            print 'Worker No. %d. is executing %s on miRNA %s.' % (idx, 
                                                                   tool['tool'],
                                                                   mir_obj['identifier'])
        pre_time = datetime.datetime.now()
        
        miRmap_out = ''
        # For each list of targets
        for j in xrange(0, len(targets_list)):
            target = targets_list[j]
            miRmap_out += '>%s %s\n\n' % (mir_obj['identifier'], target['identifier'])
            miRmap_out += self.mirmapRuns(mir_obj['seq'], target['seq'])
            miRmap_out += '\n\n' 
        
        with open(output_path, 'w') as f:
            f.write(miRmap_out)
        
        next_time = datetime.datetime.now()
        
        delta_time = next_time - pre_time
        running_time = int(delta_time.total_seconds() * 10**6)
        
        if self.debug:
            if os.path.exists(output_path):
                print '\t*** %s has been executed on miRNA %s (Worker No. %d, Time \
                (microseconds): %d)' % (tool['tool'], 
                                        mir_obj['identifier'], 
                                        idx,
                                        running_time)
            else:
                print '\t*** An error occurred during %s execution \
                on miRNA %s (Worker No. %d)' %  (tool['tool'], mir_obj['identifier'], idx) 
                
        queue.put(({ 
            'tool': tool['tool'], 
            'miRNA': mir_obj['identifier'], 
            'time': running_time,
            'error': '' 
        }, idx))
            
        return
    
    
    def getTargetingToolsList(self):
        
        tools = {
            'miRmap': { 
                'tool': 'miRmap', 
                'cmd': 'mirmap_runs.py -m %s -n %s -t %s -i %s',
                'output_filename': 'output.txt',
                'output_splitted_filename': 'output_part_%d.txt',
                'additional_opts': ''
            }, 
            'PITA': { 
                'tool': 'PITA', 
                'cmd': 'cd %s && pita_prediction.pl -mir %s -utr %s %s -prefix %s/%s >/dev/null 2>&1',
                'output_filename': 'output',
                'output_splitted_filename': 'output_part_%d',
                'additional_opts': '-l "7-8" -gpx -gu "7;0;8;0" -m "7;0;8;0"'
            }, 
            'TargetScan': { 
                'tool': 'TargetScan', 
                'cmd': 'targetscan_70.pl %s %s %s',
                'cmd2': 'targetscan_70_BL_PCT.pl %s %s %s > %s 2>/dev/null',
                'output_filename': 'output_sort.txt',
                'output_filename2': 'output.txt',
                'output_splitted_filename': 'output_sort_part_%d.txt',
                'output_splitted_filename2': 'output_part_%d.txt',
                'additional_opts': ''
            },
            'miRanda': { 
                'tool': 'miRanda', 
                'cmd': 'miranda %s %s %s -out %s/%s -en -20',
                'output_filename': 'output.txt',
                'output_splitted_filename': 'output_part_%d.txt',
                'additional_opts': '-quiet'
            },
            'RNAhybrid': { 
                'tool': 'RNAhybrid', 
                'cmd': 'RNAhybrid -c -e -20 -s 3utr_human -q %s -t %s %s > %s/%s',
                'output_filename': 'rnahybrid_output.txt',
                'output_splitted_filename': 'output_part_%d.txt',
                'additional_opts': '-m %d -n %d'            
            }
        }
        
        return tools
    
    
    def waitAndGetWorkerResults(self, jobs, queue, results_list):
        
        # Collects all results into in a list                    
        for j in xrange(0, len(jobs)):
            
            # Get worker result
            res = queue.get()
            if res:
                tool_name = res[0]['tool']
                running_time = res[0]['time']
                error = res[0]['error']
                results_list[tool_name].append({
                    'time': running_time,
                    'error': error
                })
    
    
    def waitAndGetWorkerFinalResults(self, jobs, queue, results_list):
        
        # Collects all results into in a list                    
        for j in xrange(0, len(jobs)):
            
            # Get worker result
            res = queue.get()
            if res:
                tool_name = res[0]['tool']
                results = res[0]['results']
                results_list[tool_name].append(results)
                
        # Wait for all worker processes to finish
        for p in jobs:
            p.join()    
            
    def executeTargetScan(self, tool, mirna_obj, targets_obj, jobs, queue, 
                          exec_type):
               
        if exec_type == 'multiple':
            # TargetScan
            tmp_targets_split_list = [[] for i in xrange(0, self.no_cores)]
            
            for i in xrange(0, self.no_cores):
                for j in xrange(0, len(targets_obj.targets_targetscan_split_details)):
                    idx = (j*self.no_cores) + i
                    if idx < len(targets_obj.targets_targetscan_split_details):
                        tmp_targets_split_list[i].append(targets_obj.targets_targetscan_split_details[idx])
                        
            miRNA_tool_dir_path = '%s/%s' % (mirna_obj['output_path'], 
                                             tool['tool']) 
            
            # miRNA FASTA file path
            miRNA_f_path = '%s/mirna_targetscan.txt' % mirna_obj['output_path']
            
            # Create the output tool directory
            self.createDirectory(miRNA_tool_dir_path)
            
            if os.path.exists(miRNA_tool_dir_path) and \
               os.path.exists(miRNA_f_path):            

                # For each pool of splitted files
                for s in xrange(0, len(tmp_targets_split_list)):              
                    
                    p = mp.Process(target = self.TargetScanWorker, 
                                   args=(tool, 
                                         s, 
                                         mirna_obj,
                                         tmp_targets_split_list[s],
                                         queue,
                                         exec_type))
                    jobs.append(p)
                    p.start()
                    
                
    def executeMiRmap(self, tool, mirna_obj, targets_obj, jobs, queue, exec_type):
        if exec_type == 'single':
            # Output file path
            output_f_path = tool['output_filename']
            
            # miRNA FASTA file path
            miRNA_f_path = '%s/mirna.fa' % mirna_obj['output_path']
            
            miRNA_tool_dir_path = '%s/%s' % (mirna_obj['output_path'], 
                                            tool['tool'])    
            
            miRmap_output = '%s/%s' % (miRNA_tool_dir_path, output_f_path)
            
            # Create the output tool directory
            self.createDirectory(miRNA_tool_dir_path)                                        
            
            if os.path.exists(miRNA_tool_dir_path) and \
               os.path.exists(miRNA_f_path):
                p = mp.Process(target = self.miRmapWorker, 
                               args=(0, 
                                     tool,
                                     mirna_obj,
                                     targets_obj.targets_list,
                                     miRmap_output,
                                     queue))
                jobs.append(p)
                p.start()
        else:
            # For each splitted file
            for s in xrange(0, len(targets_obj.targets_split_list)):
                
                # Output file path
                output_f_path = tool['output_splitted_filename'] % s
                
                # miRNA FASTA file path
                miRNA_f_path = '%s/mirna.fa' % mirna_obj['output_path']
                
                miRNA_tool_dir_path = '%s/%s' % (mirna_obj['output_path'], 
                                                tool['tool'])    
                
                miRmap_output = '%s/%s' % (miRNA_tool_dir_path, output_f_path)
                
                # Create the output tool directory
                self.createDirectory(miRNA_tool_dir_path)                                        
                
                if os.path.exists(miRNA_tool_dir_path) and \
                   os.path.exists(miRNA_f_path):
                    p = mp.Process(target = self.miRmapWorker, 
                                   args=(s, 
                                         tool,
                                         mirna_obj,
                                         targets_obj.targets_split_list[s],
                                         miRmap_output,
                                         queue))
                    jobs.append(p)
                    p.start()                                        
                                
    
    def executeRnahybridPita64bitMirandaTools(self, tool, mirna_obj, 
                                              targets_obj, jobs, queue, 
                                              exec_type):
        if exec_type == 'single':
            # Output file path
            output_f_path = tool['output_filename']
            
            # miRNA FASTA file path
            miRNA_f_path = '%s/mirna.fa' % mirna_obj['output_path']
            
            miRNA_tool_dir_path = '%s/%s' % (mirna_obj['output_path'], 
	                                     tool['tool'])      
            
            # Create the output tool directory
            self.createDirectory(miRNA_tool_dir_path)                            
            
            if os.path.isdir(miRNA_tool_dir_path) and \
               os.path.isfile(miRNA_f_path):
                additional_opt = tool['additional_opts']

		command = ''
		
		# Only for PITA
		if tool['tool'] == 'PITA':
		    # Create the i-th output directory
		    ith_dir_path = '%s/%s' % (miRNA_tool_dir_path, 
		                              tool['output_splitted_filename'] % 0)
		    self.createDirectory(ith_dir_path) 	
		    output_f_path = tool['output_splitted_filename'] % 0
		    
		    command = tool['cmd'] % (ith_dir_path,
		                             miRNA_f_path, 
			                     targets_obj.path,
			                     additional_opt,
			                     miRNA_tool_dir_path,
			                     output_f_path)
		else:
		    # Only for RNAhybrid
		    if tool['tool'] == 'RNAhybrid':
			additional_opt = tool['additional_opts'] \
			    % (targets_obj.longest_target_seq, len(mirna_obj['seq']))
		    
		    command = tool['cmd'] % (miRNA_f_path, 
			                     targets_obj.path,
			                     additional_opt,
			                     miRNA_tool_dir_path,
			                     output_f_path)				    
		    
                p = mp.Process(target = self.worker, 
                               args=(command, 
                                     0, 
                                     tool,
                                     mirna_obj,
                                     queue))
                jobs.append(p)
                p.start()   
        else:
            # For each splitted file
            for s in xrange(0, len(targets_obj.target_split_details)):
                
                # Output file path
                output_f_path = tool['output_splitted_filename'] % s
                                
                # miRNA FASTA file path
                miRNA_f_path = '%s/mirna.fa' % mirna_obj['output_path']
                
                miRNA_tool_dir_path = '%s/%s' % (mirna_obj['output_path'], 
		                                 tool['tool'])    
                # Create the output tool directory
                self.createDirectory(miRNA_tool_dir_path)
                
                if os.path.exists(miRNA_tool_dir_path) and \
                   os.path.exists(miRNA_f_path):
                    additional_opt = tool['additional_opts']

		    command = ''
		    
		    # Only for PITA
		    if tool['tool'] == 'PITA':
		        # Create the i-th output directory
		        ith_dir_path = '%s/%s' % (miRNA_tool_dir_path, 
		                                  tool['output_splitted_filename'] %s)
		        self.createDirectory(ith_dir_path) 	
		    
		        command = tool['cmd'] % (ith_dir_path,
		                                 miRNA_f_path, 
			                         targets_obj.target_split_details[s]['path'],
			                         additional_opt,
			                         miRNA_tool_dir_path,
			                         output_f_path)	
		    else:
			# Only for RNAhybrid
			if tool['tool'] == 'RNAhybrid':
			    additional_opt = tool['additional_opts'] \
				% (targets_obj.longest_target_seq, len(mirna_obj['seq']))
			
			command = tool['cmd'] % (miRNA_f_path, 
			                         targets_obj.target_split_details[s]['path'],
			                         additional_opt,
			                         miRNA_tool_dir_path,
			                         output_f_path)			
                
                    p = mp.Process(target = self.worker, 
                                   args=(command, 
                                         s, 
                                         tool,
                                         mirna_obj,
                                         queue))
                    jobs.append(p)
                    p.start()
    
    
    def parseMirandaResults(self, mirna_obj, targets_obj, idx, output_f_path, 
                            tool, queue):
        results = []
        if os.path.exists(output_f_path):
            with open(output_f_path, 'r') as f:
                lines = f.readlines()
                for i in xrange(0, len(lines)):
                    # Query
                    matchObj = re.match(r'^\s+Query:\s+3\'\s+([^\s]+)\s+5\'$', lines[i], re.M|re.I)
                    if matchObj:
                        i += 1
                        # Pairing
                        if i < len(lines):
                            # Get the seed region 2-7
                            seed_region = lines[i][-9:-2]
                            if len(seed_region.replace('|', '')) == 0:
                                i += 1
                                if i < len(lines):
                                    matchObj1 = re.match(r'^\s+Ref:\s+5\'\s+([^\s]+)\s+3\'$', lines[i], re.M|re.I)
                                    i += 5
                                    # miRNA - Target
                                    if i < len(lines):
                                        matchObj2 = re.match(r'^>([^\s]+)\s*.*([NMRX]{2}_[0-9_\.]+)\s+.*$', lines[i], re.M|re.I)
                                        if matchObj2:
                                            tar = matchObj2.group(2)
                                            if tar not in results:
                                                results.append(tar)
        queue.put(({ 
            'tool': tool['tool'], 
            'miRNA': mirna_obj['identifier'], 
            'results': results 
        }, idx))
            
        return
    
    
    def parseMirmapResults(self, mirna_obj, targets_obj, idx, output_f_path, tool, queue):
        results = []
        if os.path.exists(output_f_path):
            with open(output_f_path, 'r') as f:
                lines = f.readlines()
                for i in xrange(0, len(lines)):
                    # miRNA - Target
                    matchObj = re.match(r'^>([^\s]+)\s*.*([NMRX]{2}_[0-9_\.]+)$', lines[i], re.M|re.I)
                    if matchObj:
                        tar = matchObj.group(2)
                        i += 2
                        if i < len(lines):
                            matchObj = re.match(r'^[0-9]+.*$', lines[i], re.M|re.I)
                            if matchObj:
                                if tar not in results:
                                    results.append(tar)
        queue.put(({ 
            'tool': tool['tool'], 
            'miRNA': mirna_obj['identifier'], 
            'results': results 
        }, idx))
            
        return    
    
    
    def parseRnahybridResults(self, mirna_obj, targets_obj, idx, output_f_path, 
                              tool, queue):
        results = []
        if os.path.exists(output_f_path):
            with open(output_f_path, 'r') as f:        
                handler = csv.reader(f, delimiter=':')
                for line in handler:  
                    if len(line) == 11:
                        # Target
                        matchObj = re.match(r'^.*([NM]{2}_[0-9_\.]+)$', line[0], re.M|re.I)
                        if matchObj:
                            tar = matchObj.group(1) 
                            # Get the seed region 2-7
                            target_seq = line[8][-8:-1]
                            mirna_seq = line[9][-8:-1]
                            if not re.search(r'\s', mirna_seq):
                                seq_check = False
                                for i in xrange(0, len(mirna_seq)):
                                    # miRNA nucleotide
                                    mir_nt = mirna_seq[i]
                                    # Target nucleotide
                                    tar_nt = target_seq[i]
                                    if (mir_nt.upper() == 'U' and tar_nt.upper() == 'G') or \
                                       (mir_nt.upper() == 'G' and tar_nt.upper() == 'U'):
                                        seq_check = True
                                        break
                                    elif (mir_nt.upper() == 'T' and tar_nt.upper() == 'G') or \
                                       (mir_nt.upper() == 'G' and tar_nt.upper() == 'T'):
                                        seq_check = True
                                        break
                                if not seq_check and tar not in results:
                                    results.append(tar)                                 
        queue.put(({ 
            'tool': tool['tool'], 
            'miRNA': mirna_obj['identifier'], 
            'results': results 
        }, idx))
            
        return      
    
    
    def parsePita64bitResults(self, mirna_obj, targets_obj, idx, output_f_path, 
                              tool, queue):
        results = []
        if os.path.exists(output_f_path):
            with open(output_f_path, 'r') as f:        
                handler = csv.reader(f, delimiter='\t')
                for line in handler: 
                    if len(line) == 13:
                        # Target
                        matchObj = re.match(r'^.*([NMRX]{2}_[0-9_\.]+).*$', line[0], re.M|re.I)
                        if matchObj:
                            tar = matchObj.group(1) 
                            ddG = float(line[12])
                            if ddG <= -10.0:
                                if tar not in results:
                                    results.append(tar)  
        queue.put(({ 
            'tool': tool['tool'], 
            'miRNA': mirna_obj['identifier'], 
            'results': results 
        }, idx))
            
        return 
    
    
    def parseTargetscanResults(self, mirna_obj, targets_obj, idx, output_f_path, 
                               tool, queue):
        results = []
        if os.path.exists(output_f_path):
            with open(output_f_path, 'r') as f:        
                handler = list(csv.reader(f, delimiter='\t'))
                for i in xrange(1, len(handler)): 
                    line = handler[i]
                    if len(line) == 14:
                        if line[8] != '6mer':
                            # Target
                            enst_id = re.sub(r'\.[0-9]+$', '', line[0])
                            if enst_id in targets_obj.targets_enst_to_refseq_map:
                                tar = targets_obj.targets_enst_to_refseq_map[enst_id]                                
				if tar != '' and tar not in results:
				    results.append(tar)  
        queue.put(({ 
            'tool': tool['tool'], 
            'miRNA': mirna_obj['identifier'], 
            'results': results 
        }, idx))
        
        return 
    
                    
    def executeOutputParsing(self, tool, mirna_obj, targets_obj, jobs, queue, 
                             exec_type, worker):   
        if exec_type == 'single':
            # Output file path
            output_f_path = '%s/%s/%s' % (mirna_obj['output_path'], 
                                          tool['tool'],
                                          tool['output_filename'])  
		
            if tool['tool'] == 'PITA':
                # Output file path
                output_f_path = '%s/%s/%s_pita_results.tab' % (mirna_obj['output_path'], 
		                                               tool['tool'],
		                                               tool['output_splitted_filename'] % 0)
                
            p = mp.Process(target = worker, 
                           args=(mirna_obj,
                                 targets_obj,
                                 0, 
                                 output_f_path,
                                 tool,
                                 queue))
            jobs.append(p)
            p.start()  
        else:
            targets_list = []
            if tool['tool'] == 'TargetScan':
                targets_list = targets_obj.targets_targetscan_split_details
            else:
                targets_list = targets_obj.target_split_details
            
            # For each splitted file
            for s in xrange(0, len(targets_list)):
                # Output file path
                output_f_path = '%s/%s/%s' % (mirna_obj['output_path'], 
                                              tool['tool'],
                                              tool['output_splitted_filename'] % s)
                
                if tool['tool'] == 'PITA':
                    # Output file path
		    output_f_path = '%s/%s/%s_pita_results.tab' % (mirna_obj['output_path'], 
		                                                   tool['tool'],
		                                                   tool['output_splitted_filename'] % s)  
                    
                if tool['tool'] == 'TargetScan':
                    # Output file path
                    output_f_path = '%s/%s/%s' % (mirna_obj['output_path'], 
		                                  tool['tool'],
		                                  tool['output_splitted_filename2'] % s)
                                    
                p = mp.Process(target = worker, 
                               args=(mirna_obj, 
                                     targets_obj,
                                     s, 
                                     output_f_path,
                                     tool,
                                     queue))               
                jobs.append(p)
                p.start()                          
    
    
    def executeOutputMerging(self, tool, mirna_obj, targets_obj, jobs, queue, 
                             exec_type, worker):   
        if exec_type == 'multiple':
            # Output file path
            output_f_path = '%s/%s/%s' % (mirna_obj['output_path'], 
                                          tool['tool'],
                                          tool['output_filename'])  
		
            if tool['tool'] == 'PITA':
                # Output file path
                output_f_path = '%s/%s/%s_pita_results.tab' % (mirna_obj['output_path'], 
		                                               tool['tool'],
		                                               tool['output_splitted_filename'] % 0)
                
            if tool['tool'] == 'TargetScan':
                # Output file path
                output_f_path = '%s/%s/%s' % (mirna_obj['output_path'],
		                              tool['tool'],
		                              tool['output_filename2']) 
                
            p = mp.Process(target = worker, 
                           args=(mirna_obj,
                                 targets_obj,
                                 0, 
                                 output_f_path,
                                 tool,
                                 queue))
            jobs.append(p)
            p.start()  
        else:
            targets_list = []
            targets_list = targets_obj.target_split_details
            
            # For each splitted file
            for s in xrange(0, len(targets_list)):
                # Output file path
                output_f_path = '%s/%s/%s' % (mirna_obj['output_path'], 
                                              tool['tool'],
                                              tool['output_splitted_filename'] % s)
                
                if tool['tool'] == 'PITA':
                    # Output file path
		    output_f_path = '%s/%s/%s_pita_results.tab' % (mirna_obj['output_path'], 
		                                                   tool['tool'],
		                                                   tool['output_splitted_filename'] % s)  
               
                p = mp.Process(target = worker, 
                               args=(mirna_obj, 
                                     targets_obj,
                                     s, 
                                     output_f_path,
                                     tool,
                                     queue))               
                jobs.append(p)
                p.start()                   
    
    def executeTargetingAnalysis(self, tools_list, selected_tools_list, 
                                 no_runs, miRNAs_obj, targets_obj, queue, 
                                 out_queue):
        results_list = {}
        final_results_list = {} 
	
	# Contain WT, NT substitution and isomiR miRNAs
	all_mirnas_list = []
	out_merge_list = []
	
	for mir in miRNAs_obj.mirnas_list:
	   
	    all_mirnas_list.append({
	        'identifier': mir['identifier'],
	        'output_path': mir['output_path'],
	        'length': mir['length'],
	        'seq': mir['seq']
	    })
	    
	    for m in mir['variations']['substitution']:
		
		all_mirnas_list.append({
		    'identifier': m['identifier'],
		    'output_path': m['output_path'],
		    'length': m['length'],
		    'seq': m['seq']
		})	
		
	    for m in mir['variations']['isomir']:
		
		all_mirnas_list.append({
		    'identifier': m['identifier'],
		    'output_path': m['output_path'],
		    'length': m['length'],
		    'seq': m['seq']
		})
	
        # For each miRNA
        for mir in all_mirnas_list:
            results_list[mir['identifier']] = {}
            
            # For each tool
            for tl in selected_tools_list:
                results_list[mir['identifier']][tl] = []
                jobs = []
                if tl in tools_list:
                    # Get tool details
                    s_tool = tools_list[tl]
                    if tl in ['miRanda', 'RNAhybrid', 'PITA']:
                        # ======================================================
                        # miRanda, RNAhybrid, PITA
                        # ======================================================                                       
                        # Multiple execution
                        if len(targets_obj.target_split_details) > 0:
                            self.executeRnahybridPita64bitMirandaTools(
                                s_tool, mir, targets_obj, jobs, queue, 'multiple')
                        else:
                            # Single execution
                            self.executeRnahybridPita64bitMirandaTools(
                                s_tool, mir, targets_obj, jobs, queue, 'single')
                    elif tl == 'miRmap':
                        # ======================================================
                        # miRmap
                        # ======================================================                   
                        # Multiple execution
                        if len(targets_obj.target_split_details) > 0:
                            self.executeMiRmap(s_tool, mir, targets_obj, 
                                               jobs, queue, 'multiple')
                        else:
                            # Single execution
                            self.executeMiRmap(s_tool, mir, targets_obj, 
                                               jobs, queue, 'single')
                    elif tl == 'TargetScan':
                        # ======================================================
                        # TargetScan
                        # ======================================================                        
                        # Multiple execution  
                        if len(targets_obj.targets_targetscan_split_details) > 1:
                            self.executeTargetScan(s_tool, mir, targets_obj, 
                                                   jobs, queue, 'multiple')
                        
                self.waitAndGetWorkerResults(jobs, queue, results_list[mir['identifier']]) 
               
            #===================================================================
            #===================================================================
            # Tools output analysis
            #===================================================================
            #===================================================================
            # Output file path
            final_output_path = '%s/output.txt' % mir['output_path'] 
            # Output file path (mirwalk format)
            mirwalk_final_output_path = '%s/output.tsv' % mir['output_path']             
            
            final_results_list[mir['identifier']] = {}
            
            # Store all unique targets among selected tools 
            all_targets = []
            
            # For each tool 
            for tl in selected_tools_list:
                final_results_list[mir['identifier']][tl] = []
                jobs = []
                if tl in tools_list:
                    # Get tool details
                    s_tool = tools_list[tl] 
                    
                    # ==========================================================
                    # miRanda
                    # ==========================================================  
                     # Check if the tool is available
                    if tl == 'miRanda':
                        # Multiple execution
                        if len(targets_obj.target_split_details) > 0:
                            self.executeOutputParsing(s_tool, mir, targets_obj, 
                                                      jobs, out_queue, 'multiple',
                                                      self.parseMirandaResults)
                        else:
                            # Single execution
                            self.executeOutputParsing(s_tool, mir, targets_obj, 
                                                      jobs, out_queue, 'single',
                                                      self.parseMirandaResults)
                    elif tl == 'miRmap':
                        # ======================================================
                        # miRmap
                        # ====================================================== 
                        # Multiple execution
                        if len(targets_obj.target_split_details) > 0:
                            self.executeOutputParsing(s_tool, mir, targets_obj, 
                                                      jobs, out_queue, 'multiple',
                                                      self.parseMirmapResults)
                        else:
                            # Single execution
                            self.executeOutputParsing(s_tool, mir, targets_obj,
                                                      jobs, out_queue, 'single',
                                                      self.parseMirmapResults)
                    elif tl == 'RNAhybrid':
                        # ======================================================
                        # RNAhybrid
                        # ====================================================== 
                        # Multiple execution
                        if len(targets_obj.target_split_details) > 0:
                            self.executeOutputParsing(s_tool, mir, targets_obj, 
                                                      jobs, out_queue, 'multiple',
                                                      self.parseRnahybridResults)
                        else:
                            # Single execution
                            self.executeOutputParsing(s_tool, mir, targets_obj,
                                                      jobs, out_queue, 'single',
                                                      self.parseRnahybridResults)
                    elif tl == 'PITA':
                        # ======================================================
                        # PITA
                        # ====================================================== 
                        # Multiple execution
                        if len(targets_obj.target_split_details) > 0:
                            self.executeOutputParsing(s_tool, mir, targets_obj, 
                                                      jobs, out_queue, 'multiple',
                                                      self.parsePita64bitResults)
                        else:
                            # Single execution
                            self.executeOutputParsing(s_tool, mir, targets_obj,
                                                      jobs, out_queue, 'single',
                                                      self.parsePita64bitResults) 
                    elif tl == 'TargetScan':
                        # ======================================================
                        # TargetScan
                        # ====================================================== 
                        # Multiple execution
                        if len(targets_obj.targets_targetscan_split_details) > 1:
                            self.executeOutputParsing(s_tool, mir, targets_obj, 
                                                      jobs, out_queue, 'multiple',
                                                      self.parseTargetscanResults)
                       
                    self.waitAndGetWorkerFinalResults(
                        jobs, out_queue, final_results_list[mir['identifier']])                    
                       
                    # Merge all found targets list according to a specific miRNA and tool
                    merged_targets = list(set().union(*final_results_list[mir['identifier']][tl]))
                    merged_targets.sort()
                    final_results_list[mir['identifier']][tl] = list(set(merged_targets))
                    
		    #all_targets = all_targets + [ x for x in merged_targets if x not in all_targets ]
		    all_targets += merged_targets

	    # Remove all duplicates
	    all_targets = list(set(all_targets))
            all_targets.sort()
	    
	    # matrix stores targets according to the number of tools that identify each target
            matrix = {
	        1: [],
	        2: [],
	        3: [],
	        4: [],
	        5: []
	    }
            
            with open(final_output_path, 'w') as f:
                f.write('Tool\tNo. targets\tTargets\n')
                for tool,results in final_results_list[mir['identifier']].iteritems():
                    f.write('%s\t%d\t%s\n' % (tool, len(results), ';'.join(results)))
                    
            # Generate the structure for the mirwalk format       
            for target in all_targets:
                row = {
                    'target': target,
                    'tools': {},
                    'sum': 0
                }
                
                for tool,results in final_results_list[mir['identifier']].iteritems():
                    predicted = 1 if target in results else 0
                    row['sum'] += predicted
                    row['tools'][tool] = str(predicted)
                
                matrix[row['sum']].append(row)
		
            with open(mirwalk_final_output_path, 'w') as f:
		row_tools_list = []
		for i in reversed(xrange(1, 6)):
		    if len(matrix[i]) > 0:
			row_tools_list = matrix[i][0]['tools'].keys()
			row_tools_list.sort()
			f.write('Target\t%s\tSum\n' % '\t'.join(row_tools_list))
			break
		for i in reversed(xrange(1, 6)):
		    for row in matrix[i]:
			f.write('%s\t%s\t%d\n' % (row['target'], '\t'.join([row['tools'][key] for key in row_tools_list]), row['sum']))
			
	    #===================================================================
            #===================================================================
            # Tools output merging
            #===================================================================
            #===================================================================
	    # For each tool 
	    for tl in selected_tools_list:
		jobs = []
		if tl in tools_list:
		    # Get tool details
		    s_tool = tools_list[tl] 
		    n_files = 0
		    
		    if s_tool['tool'] == 'TargetScan':
			n_files = len(targets_obj.targets_targetscan_split_details)
		    else:
			n_files = len(targets_obj.target_split_details)	    
		    
		    if n_files > 0:
			out_merge_list.append({
			    'tool': s_tool['tool'],
			    'n_files': n_files,
			    'dir_path': '{0}/{1}'.format(mir['output_path'], s_tool['tool'])
			})
    
	# Available only for multiple exectutions
	if len(out_merge_list) > 0:
	    out_dirs_per_core = [ [] for i in xrange(0, self.no_cores) ]
	    
	    # Prepare all output prediction dirs to be assigned to each available core
	    for i in xrange(0, len(out_merge_list)):
		idx = int(math.fmod(i, self.no_cores))
		# Add the i-th prediction to the idx-th core
		out_dirs_per_core[idx].append(out_merge_list[i])
	    
	    jobs = []
	    
	    # For each core
	    for i in xrange(0, len(out_dirs_per_core)):
		if out_dirs_per_core[i]:
		    p = mp.Process(target = self.mergeWorker, 
			           args=(i,
			                 out_dirs_per_core[i]))               
		    jobs.append(p)
		    p.start() 
		    
	    # Wait for all worker processes to finish
	    for p in jobs:
		p.join() 	    
	  
        return results_list, final_results_list
    
    
    
