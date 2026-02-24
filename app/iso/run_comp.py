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
import isotarpkg_comp as iso
import os
import os.path
import json

if __name__ == '__main__':
    
    # Resource paths
    output_dir = '/app/iso/static/public/output'
    go_path = '/app/iso/static/public/resources/go.json' 
    job_go_comp_path = '/app/iso/static/public/output/job_go_comp.json'
    nm_to_geneid_path = '/app/iso/static/public/resources/nm_to_geneid.tab'
    obo_dag_path = '/app/iso/static/public/resources/go-basic.obo'
    gene_2_go_path = '/app/iso/static/public/resources/gene2go'     
    gene_2_pathway_path = '/app/iso/static/public/resources/genename_pathways.json' 
    pathways_path = '/app/iso/static/public/resources/pathways.json'    
    
    isotar = iso.Isotar(output_dir,
                        go_path,
                        nm_to_geneid_path,
                        obo_dag_path,
                        gene_2_go_path,
                        pathways_path,
                        gene_2_pathway_path)
    isotar.run()
 
    # Check for a previous job go comp
    if os.path.isfile(job_go_comp_path):
	with open(job_go_comp_path, 'r') as f:
	    job = json.load(f) 
	
	    if 'job_status' in job:
		job['job_status'] = 'done'
		job['job_pid'] = -1
		# Update
		with open(job_go_comp_path, 'w') as f:
		    f.write(json.dumps(job, indent=4))
    