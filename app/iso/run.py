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
import isotarpkg as iso
import os
import os.path
import json
import datetime

if __name__ == '__main__':
    
    # Resource paths
    mirna_path = '/app/iso/static/public/output/mirnas.json'
    target_path = '/app/iso/static/public/resources/3utr.fa'
    output_dir = '/app/iso/static/public/output'
    go_path = '/app/iso/static/public/resources/go.json' 
    map_path = '/app/iso/static/public/resources/targets_go_map.tab'
    pre_mature_path = '/app/iso/static/public/resources/marute_pre_mirna.json' 
    job_path = '/app/iso/static/public/output/job.json'
    job_go_path = '/app/iso/static/public/output/job_go.json'
    pre_mat_acc_path = '/app/iso/static/public/resources/mat_pre_acc.json'
    nm_to_geneid_path = '/app/iso/static/public/resources/nm_to_geneid.tab'
    obo_dag_path = '/app/iso/static/public/resources/go-basic.obo'
    gene_2_go_path = '/app/iso/static/public/resources/gene2go'     
    gene_2_pathway_path = '/app/iso/static/public/resources/genename_pathways.json' 
    pathways_path = '/app/iso/static/public/resources/pathways.json'
    
    pre_time = datetime.datetime.now()
    
    isotar = iso.Isotar(mirna_path,
                        target_path,
                        output_dir,
                        map_path,
                        go_path,
                        pre_mature_path,
                        pre_mat_acc_path,
                        nm_to_geneid_path,
                        obo_dag_path,
                        gene_2_go_path,
                        pathways_path,
                        gene_2_pathway_path)
    isotar.run()
        
    # Check for a previous job
    if os.path.isfile(job_path):
	with open(job_path, 'r') as f:
	    job = json.load(f) 
	
	    if 'job_status' in job:
		job['job_status'] = 'done'
		job['job_pid'] = -1
		# Update
		with open(job_path, 'w') as f:
		    f.write(json.dumps(job, indent=4))  
		    
    # Check for a previous job go
    if os.path.isfile(job_go_path):
	with open(job_go_path, 'r') as f:
	    job = json.load(f) 
	
	    if 'job_status' in job:
		job['job_status'] = 'done'
		job['job_pid'] = -1
		# Update
		with open(job_go_path, 'w') as f:
		    f.write(json.dumps(job, indent=4))
		    
    next_time = datetime.datetime.now()
        
    delta_time = next_time - pre_time
    running_time = int(delta_time.total_seconds() * 10**6)
    print "Exectution time:\t", running_time
    