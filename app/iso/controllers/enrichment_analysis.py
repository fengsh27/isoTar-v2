# -*- coding: UTF-8 -*-
# Import flask dependencies
from flask import (
    Blueprint, 
    render_template, 
    make_response,
    jsonify,
    request
)
from jinja2 import TemplateNotFound
import requests
import json
from iso.utilities import (
    search,
    getEditingSitesList,
    getBoxPlotData,
    getColorPicker,
    getRGBAColorPicker,
    getRepositoryVersion
)
import re
import psutil
import json
import os
import os.path
from iso.settings import AppResources as res
import iso.isotarpkg as iso
import time
import subprocess
from iso.tasks import celery_run_command
import time
from marshmallow import ValidationError
import signal
import os, errno

enrichment_analysis = Blueprint('enrichment_analysis', __name__)


def available_cpu_count():
    """ Number of available virtual or physical CPUs on this system, i.e.
    user/real as output by time(1) when called with an optimally scaling
    userspace-only program"""

    # https://github.com/giampaolo/psutil
    try:
	return psutil.cpu_count()   # psutil.NUM_CPUS on old versions
    except (ImportError, AttributeError):
	pass

    # POSIX
    try:
	res = int(os.sysconf('SC_NPROCESSORS_ONLN'))
	if res > 0:
	    return res
    except (AttributeError, ValueError):
	pass

    # Windows
    try:
	res = int(os.environ['NUMBER_OF_PROCESSORS'])
	if res > 0:
	    return res
    except (KeyError, ValueError):
	pass

    return 1


def execute_command_line(command):
    retcode = 0
    try:
	retcode = subprocess.check_output(command, shell = True)
    except OSError as e:
	print >>sys.stderr, "Execution failed:", e   
	
    return retcode 

def remove_file(path):
    if os.path.isfile(path):
	try:
	    os.remove(path)
	except OSError as e:
	    if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
		pass
	
	
@enrichment_analysis.route('/enrichment_analysis')
def index_enrichment_analysis():
    
    # Get the number of CPU cores
    #no_system_cores = psutil.cpu_count()
    no_system_cores = available_cpu_count()
    
    """Response"""    
    try:
        return render_template('enrichment_analysis/enrichment.html',
	                       sys_cores = no_system_cores)
    except TemplateNotFound:
        return make_response('Not Found'), 404
    

@enrichment_analysis.route('/exec_functional_enrichment_analysis', methods=['POST'])
def exec_functional_enrichment_analysis():
    
    if request.method != 'POST':
        return jsonify({ 'error_msg': 'No data available' }), 405
    
    if os.path.isfile(res.JOB_PATH) or \
       os.path.isfile(res.JOB_ENRICH_PATH) or \
       os.path.isfile(res.JOB_ENRICH_COMP_PATH):
	return jsonify({ 'error_msg': 'There is already a running job.' }), 405    
   
    if not os.path.isfile(res.PRE_JSON_PATH):
	return jsonify({ 'error_msg': 'No data available' }), 405  
    
    if not os.path.isfile(res.TASK_PATH):
	return jsonify({ 'error_msg': 'No data available' }), 405     
    
    cores = request.form.getlist('cores')[0]
    pvalue = request.form.getlist('pvalue')[0]
    pvalue_type = request.form.getlist('pvalue_type')[0]
    enrich_opts = request.form.getlist('enrich_opts[]')
    
    # ==========================================================================
    # Input filtering 
    # ==========================================================================
    cores_f = int(re.sub(r'[^0-9]+', '', cores))
    pvalue_f = float(re.sub(r'[^0-9\.]+', '', pvalue))
    pvalue_type_f = re.sub(r'[^a-z]+', '', pvalue_type)
    consensus = 1
    enrich_opts_f = []
    
    for el in enrich_opts:
        f_el = re.sub(r'[^a-z]', '', el)
        if f_el in ['go', 'kegg', 'reactome']:
            enrich_opts_f.append(f_el) 
	    
    with open(res.TASK_PATH, 'r') as f:
	d = json.load(f)
	if 'consensus' in d:
	    consensus = int(d['consensus'])
	
    if not isinstance(cores_f, int) or \
       not isinstance(pvalue_f, float) or \
       not isinstance(consensus, int) or \
       pvalue_type_f not in ['raw', 'adj'] or \
       len(enrich_opts_f) == 0:
        return jsonify({ 'error_msg': 'No data available' }), 405  
        
    command = "/usr/bin/python /app/iso/run.py --cores {0} --enrich yes --cons {1} --pvalue {2} --pvalue_type {3} --enrichtype {4}".format(
        cores_f, consensus, pvalue_f, pvalue_type_f, ','.join(enrich_opts_f))
    
    try:
	task = celery_run_command.delay(command.split()) 
	task_id = task.id
	task_status = celery_run_command.AsyncResult(task_id)
	#task_state = task_status.state
	job_details = task.get()
	
	job = {
	    'task_id': task_id, 
	    'job_pid': job_details['pid'],
	    'job_status': job_details['status'],
	    'job_msg': job_details['msg'],
	    'command': command
	}
	
	with open(res.JOB_ENRICH_PATH, 'w') as f:
	    f.write(json.dumps(job, indent=4))	
	
	return jsonify(job), 202

    except ValidationError as e:
	return jsonify({ 'error_msg': e.messages }), 405
    
       
@enrichment_analysis.route('/get_functional_enrichment_analysis_content', methods=['POST'])
def get_functional_enrichment_analysis_content():
   
    if request.method != 'POST':
        return jsonify({ 'error_msg': 'No data available' }), 405
    
    analysis_type = request.form.getlist('analysis_type')[0]
    
    # ==========================================================================
    # Input filtering 
    # ==========================================================================
    analysis_type_f = re.sub(r'[^a-z]+', '', analysis_type)
    
    if analysis_type_f not in ['comparison', 'single']:
	return jsonify({ 'error_msg': 'No data available' }), 405
    
    f_path = res.JOB_ENRICH_PATH
    if analysis_type_f == 'comparison':  
	f_path = res.JOB_ENRICH_COMP_PATH
	
    if not os.path.isfile(f_path):
	return make_response('No data available'), 405
    
    """Response"""    
    try:
        return render_template('enrichment_analysis/enrich_results.html',
	                       analysis_type = analysis_type_f)
    except TemplateNotFound:
        return make_response('Not Found'), 404
    

@enrichment_analysis.route('/enrichment_analysis_status_check', methods=['POST'])
def enrichment_analysis_status_check():
    
    if request.method != 'POST':
        return jsonify({ 'error_msg': 'No data available' }), 405
    
    analysis_type = request.form.getlist('analysis_type')[0]
    
    # ==========================================================================
    # Input filtering 
    # ==========================================================================
    analysis_type_f = re.sub(r'[^a-z]+', '', analysis_type)
    
    if analysis_type_f not in ['comparison', 'single']:
	return jsonify({ 'error_msg': 'No data available' }), 405
    
    f_path = res.JOB_ENRICH_PATH
    if analysis_type_f == 'comparison':  
	f_path = res.JOB_ENRICH_COMP_PATH
	
    if not os.path.isfile(f_path):
	return jsonify({ 'error_msg': 'No data available' }), 405
    
    pid = None
    with open(f_path, 'r') as f:
	job = json.load(f)
	pid = job['job_pid']
	
    try:
	if pid == -1:
	    # Check the job
	    if os.path.isfile(f_path):
		with open(f_path, 'r') as f:
		    job = json.load(f)
		    job['job_status'] = 'done'
		    # Update
		    with open(f_path, 'w') as f:
			f.write(json.dumps(job, indent=4))	
	    return jsonify({ 'status': 'done' }), 200 
	
	# Damn zombie!
	proc = psutil.Process(pid)
	if proc.status() == psutil.STATUS_ZOMBIE:
	    # Check the job
	    if os.path.isfile(f_path):
		with open(f_path, 'r') as f:
		    job = json.load(f)
		    job['job_status'] = 'done'
		    job['job_pid'] = -1
		    # Update
		    with open(f_path, 'w') as f:
			f.write(json.dumps(job, indent=4))	
	    return jsonify({ 'status': 'done' }), 200 
	
	return jsonify({ 'status': 'running' }), 200 
    
    except psutil.NoSuchProcess:
	pass
    
    # Check the job
    if os.path.isfile(f_path):
	with open(f_path, 'r') as f:
	    job = json.load(f)
	    job['job_status'] = 'done'
	    job['job_pid'] = -1
	    # Update
	    with open(f_path, 'w') as f:
		f.write(json.dumps(job, indent=4))	
    return jsonify({ 'status': 'done' }), 200 
    
    
@enrichment_analysis.route('/enrichment_analysis_cleaner', methods=['POST'])
def enrichment_analysis_cleaner():

    if request.method != 'POST':
        return jsonify({ 'error_msg': 'No data available' }), 405
    
    analysis_type = request.form.getlist('analysis_type')[0]
    
    # ==========================================================================
    # Input filtering 
    # ==========================================================================
    analysis_type_f = re.sub(r'[^a-z]+', '', analysis_type)
    
    if analysis_type_f not in ['comparison', 'single']:
	return jsonify({ 'error_msg': 'No data available' }), 405
    
    f_path = res.JOB_ENRICH_PATH
    if analysis_type_f == 'comparison':  
	f_path = res.JOB_ENRICH_COMP_PATH
	
    # Delete te job_go.json file
    remove_file(f_path)
    
    # Remove all potential zombies... kill the zombies!
    process = subprocess.Popen(['bash','-c', '/opt/kill_zombies.sh'], stdout=subprocess.PIPE)
    output, error = process.communicate()

    """Response"""    
    try:
        return jsonify({ 'status': 'done' }), 200 
    except TemplateNotFound:
        return jsonify({ 'error_msg': 'Not Found' }), 404
    
    
@enrichment_analysis.route('/get_enrichment_analysis_results', methods=['POST'])
def get_enrichment_analysis_results():

    if request.method != 'POST':
        return jsonify({ 'error_msg': 'No data available' }), 405
    
    analysis_type = request.form.getlist('analysis_type')[0]
    
    # ==========================================================================
    # Input filtering 
    # ==========================================================================
    analysis_type_f = re.sub(r'[^a-z]+', '', analysis_type)
    
    if analysis_type_f not in ['comparison', 'single']:
	return jsonify({ 'error_msg': 'No data available' }), 405
    
    
    single_analysis_table_content = {
        'go': {
            'isomir': False,
            'nt_var': False,
            'wt': False
        },
        'kegg': {
            'isomir': False,
            'nt_var': False,
            'wt': False
        },
        'reactome': {
            'isomir': False,
            'nt_var': False,
            'wt': False
        }        
    }
    
    resources_list = []
    
    comparison_analysis_table_content = {
        'go': False,
        'kegg': False,
        'reactome': False
    }    
    
    enrich_data = {}
    file_path = 'enrichment_analysis/enrich_data.html'
    
    if analysis_type_f == 'single':
	# Check if show isomir/nt varietion/wt for go, kegg and reactome
	if os.path.isfile(res.ENRICH_JSON_PATH):
	    with open(res.ENRICH_JSON_PATH, 'r') as f:
		enrich_data = json.load(f)
		
	for tab,target in single_analysis_table_content.iteritems():
	    if tab in enrich_data:
		resources_list.append(tab)
		for mirna,details in enrich_data[tab].iteritems():
		    if len(details['variations']['isomir']) > 0:
			single_analysis_table_content[tab]['isomir'] = True
		    if len(details['variations']['substitution']) > 0:
			single_analysis_table_content[tab]['nt_var'] = True
		    
		    # Wild type only
		    if len(details['variations']['isomir']) == 0 and \
		       len(details['variations']['substitution']) == 0:
			single_analysis_table_content[tab]['wt'] = True
	
    elif analysis_type_f == 'comparison':
	file_path = 'enrichment_analysis/enrich_comp_data.html'
	if os.path.isfile(res.ENRICH_COMP_JSON_PATH):
	    with open(res.ENRICH_COMP_JSON_PATH, 'r') as f:
		enrich_data = json.load(f)
		
		for tab,value in comparison_analysis_table_content.iteritems():
		    if tab in enrich_data:
			resources_list.append(tab)
			comparison_analysis_table_content[tab] = True
			
    """Response"""    
    try:
        return render_template(file_path,
	                       single_analysis_table_content = single_analysis_table_content,
	                       comparison_analysis_table_content = comparison_analysis_table_content,
	                       enrich_data = enrich_data,
	                       analysis_type = analysis_type_f,
	                       resources_list = resources_list)
    except TemplateNotFound:
        return jsonify({ 'error_msg': 'Not Found' }), 404
    
    
@enrichment_analysis.route('/enrichment_analysis_variation_statistics', methods=['POST'])
def enrichment_analysis_variation_statistics():
    
    if request.method != 'POST':
        return jsonify({ 'error_msg': 'No data available' }), 405
    
    data = request.form.to_dict()
    analysis_type = request.form.getlist('analysis_type')[0]
    
    # ==========================================================================
    # Input filtering 
    # ==========================================================================
    analysis_type_f = re.sub(r'[^a-z]+', '', analysis_type)
    
    if analysis_type_f not in ['comparison', 'single']:
	return jsonify({ 'error_msg': 'No data available' }), 405
    
    table_data = []
    filename = ''
    
    if analysis_type_f == 'single':	
	
	if 'mirna_id' not in data or \
	   'vtype' not in data or \
	   'data_tab' not in data:
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	mirna_id = re.sub(r'[^0-9a-zA-Z\-]+', '', request.form.getlist('mirna_id')[0]) if request.form.getlist('mirna_id')[0] else '' 
	vtype = re.sub(r'[^\sa-zA-Z]+', '', request.form.getlist('vtype')[0]) if request.form.getlist('vtype')[0] else ''
	data_tab = re.sub(r'[^\sa-zA-Z\_]+', '', request.form.getlist('data_tab')[0]) if request.form.getlist('data_tab')[0] else ''
	data_tab = data_tab.replace('enrich_row_', '')
	
	v3p = re.sub(r'[^\-\+0-9]+', '', request.form.getlist('v3p')[0]) if request.form.getlist('v3p')[0] else ''
	v5p = re.sub(r'[^\-\+0-9]+', '', request.form.getlist('v5p')[0]) if request.form.getlist('v5p')[0] else ''
	
	vpos = re.sub(r'[^0-9]+', '', request.form.getlist('vpos')[0]) if request.form.getlist('vpos')[0] else ''
	vref = re.sub(r'[^AGCUT\-\+]+', '', request.form.getlist('vref')[0]) if request.form.getlist('vref')[0] else ''
	vmod = re.sub(r'[^AGCUT\-\+]+', '', request.form.getlist('vmod')[0]) if request.form.getlist('vmod')[0] else ''
	
	if mirna_id == '' or \
	   vtype == '' or \
	   vtype not in ['Wild Type', 'isomiR', 'Substitution'] or \
	   data_tab not in [
	       'cc', 'mf', 'bp', 
	       'cc_wt_o', 'mf_wt_o', 'bp_wt_o', 
	       'cc_ins', 'mf_ins', 'bp_ins', 
	       'cc_var_o', 'mf_var_o', 'bp_var_o'
	    ]:
	    return jsonify({ 'error_msg': 'No data available' }), 405    
	
	enrich_data = {}
	if not os.path.isfile(res.ENRICH_JSON_PATH):
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	with open(res.ENRICH_JSON_PATH, 'r') as f:
	    enrich_data = json.load(f)	
	    
	if not enrich_data:
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	if mirna_id not in enrich_data['go']:
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	stats = enrich_data['go'][mirna_id]
	
	# Check the selected GO term category and class
	data_class = ''
	data_cat = ''    
	data_label = data_tab.split('_')
	
	if len(data_label) == 1:
	    data_class = data_label[0]
	elif len(data_label) == 2:
	    data_class = data_label[0]
	    data_cat = 'intersection'
	elif len(data_label) == 3:
	    data_class = data_label[0]
	    data_cat = 'wt_only'  
	    if data_label[1] == 'var':
		data_cat = 'variation_only'
		
	if data_class != '':
	    if vtype == 'Substitution':
		if int(vpos) == -1 or vref == '' or vmod == '':
		    return jsonify({ 'error_msg': 'No data available' }), 405
		for s in stats['variations']['substitution']:
		    if s['position'] == int(vpos) and \
		       s['ref'] == vref and \
		       s['mod'] == vmod:			
			if data_cat != '':
			    if data_cat in s and data_class in s[data_cat]:
				table_data = s[data_cat][data_class]
				filename = '{0}_nt_var_{1}_{2}_{3}_{4}_go_{5}'.format(mirna_id, 
				                                                      s['position'],
				                                                      s['ref'],
				                                                      s['mod'],
				                                                      data_cat,
				                                                      data_class)				
			elif data_class in s:
			    table_data = s[data_class]
			    filename = '{0}_nt_var_{1}_{2}_{3}_go_{4}'.format(mirna_id, 
			                                                      s['position'],
			                                                      s['ref'],
			                                                      s['mod'],
			                                                      data_class)			    
			break
			
	    elif vtype == 'isomiR':
		if v3p == '' or v5p == '':
		    return jsonify({ 'error_msg': 'No data available' }), 405
		for i in stats['variations']['isomir']:
		    if i['5p'] == v5p and i['3p'] == v3p:			
			if data_cat != '':
			    if data_cat in i and data_class in i[data_cat]:
				table_data = i[data_cat][data_class]
				filename = '{0}_isomir_{1}_{2}_{3}_go_{4}'.format(mirna_id, 
				                                                  i['5p'],
				                                                  i['3p'],
				                                                  data_cat,
				                                                  data_class)				
			elif data_class in i:
			    table_data = i[data_class]
			    filename = '{0}_isomir_{1}_{2}_go_{3}'.format(mirna_id, 
			                                                  i['5p'],
			                                                  i['3p'],
			                                                  data_class)			    
			break
			    
	    elif vtype == 'Wild Type':
		if data_class in stats:
		    filename = '{0}_wild_type_go_{1}'.format(mirna_id, data_class)
		    table_data = stats[data_class]
    
    elif analysis_type_f == 'comparison':
	if 'list_id' not in data or \
	   'data_tab' not in data:
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	list_id = re.sub(r'[^0-9a-zA-Z\-]+', '', request.form.getlist('list_id')[0]) if request.form.getlist('list_id')[0] else '' 
	data_tab = re.sub(r'[^\sa-zA-Z\_]+', '', request.form.getlist('data_tab')[0]) if request.form.getlist('data_tab')[0] else ''
	data_tab = data_tab.replace('enrich_row_', '')
	
	if list_id == '' or \
	   list_id not in ['list1', 'list2'] or \
	   data_tab not in [
	       'cc', 'mf', 'bp', 
	       'cc_wt_o', 'mf_wt_o', 'bp_wt_o', 
	       'cc_ins', 'mf_ins', 'bp_ins', 
	       'cc_var_o', 'mf_var_o', 'bp_var_o'
	    ]:
	    return jsonify({ 'error_msg': 'No data available' }), 405    
	
	enrich_data = {}
	if not os.path.isfile(res.ENRICH_COMP_JSON_PATH):
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	with open(res.ENRICH_COMP_JSON_PATH, 'r') as f:
	    enrich_data = json.load(f)	
	    
	if not enrich_data:
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	if list_id not in enrich_data['go']:
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	stats = enrich_data['go'][list_id]	
	
	# Check the selected GO term category and class
	data_class = ''
	data_cat = ''    
	data_label = data_tab.split('_')
	
	if len(data_label) == 1:
	    # E.g.: bp
	    data_class = data_label[0] # bp, cc, mf
	elif len(data_label) == 2:
	    data_class = data_label[0] # bp, cc, mf
	    data_cat = 'intersection'
	    data_cat_f = 'intersection'
	elif len(data_label) == 3:
	    # E.g.: bp_wt_o
	    data_class = data_label[0] # bp, cc, mf
	    data_cat = 'wt_only'  
	    data_cat_f = 'list1_only'
	    if data_label[1] == 'var':
		data_cat = 'variation_only'
		data_cat_f = 'list2_only'
		
	if data_class != '':
	    if data_cat != '':
		if data_cat in stats and data_class in stats[data_cat]:
		    table_data = stats[data_cat][data_class]
		    filename = '{0}_{1}_go_{2}'.format(stats['entity'], 
		                                       data_cat_f,
		                                       data_class)
	    elif data_class in stats:
		table_data = stats[data_class]	
		filename = '{0}_go_{1}'.format(stats['entity'],
		                               data_tab)		
    
    """Response"""    
    try:
        return render_template('enrichment_analysis/go_statistics.html',
	                       table_data = table_data,
	                       filename = filename)
    except TemplateNotFound:
        return jsonify({ 'error_msg': 'Not Found' }), 404
    
    
@enrichment_analysis.route('/enrichment_analysis_variation_histograms', methods=['POST'])
def enrichment_analysis_variation_histograms():
    
    data = request.form.to_dict()
    analysis_type = request.form.getlist('analysis_type')[0]
    
    # ==========================================================================
    # Input filtering 
    # ==========================================================================
    analysis_type_f = re.sub(r'[^a-z]+', '', analysis_type)
    
    if analysis_type_f not in ['comparison', 'single']:
	return jsonify({ 'error_msg': 'No data available' }), 405
    
    hist_data = []
    target = ''
    pdf_filename = ''
    entity_label1 = ''
    entity_label2 = ''    
    
    section_title = ''
    
    if analysis_type_f == 'single':
	
	entity_label1 = 'Wild Type'
	entity_label2 = 'Variation'	
    
	if 'mirna_id' not in data or \
	   'vtype' not in data or \
	   'resource' not in data:
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	mirna_id = re.sub(r'[^0-9a-zA-Z\-]+', '', request.form.getlist('mirna_id')[0]) if request.form.getlist('mirna_id')[0] else '' 
	vtype = re.sub(r'[^\sa-zA-Z]+', '', request.form.getlist('vtype')[0]) if request.form.getlist('vtype')[0] else ''
	resource = re.sub(r'[^a-zA-Z]+', '', request.form.getlist('resource')[0]) if request.form.getlist('resource')[0] else ''
	
	if mirna_id == '' or \
	   vtype == '' or \
	   vtype not in ['Wild Type', 'isomiR', 'Substitution'] or \
	   resource not in ['go']:
	    return jsonify({ 'error_msg': 'No data available' }), 405    
	
	
	section_title = '{0}: {1}'.format(vtype, mirna_id)
	enrich_data = {}
	if not os.path.isfile(res.ENRICH_JSON_PATH):
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	with open(res.ENRICH_JSON_PATH, 'r') as f:
	    enrich_data = json.load(f)	
 
	if not enrich_data:
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	# GO term details
	data_details = {}
	if not os.path.isfile(res.GO_PATH):
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	with open(res.GO_PATH, 'r') as f:
	    data_details = json.load(f)
	
	if mirna_id not in enrich_data[resource]:
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	stats = enrich_data[resource][mirna_id]
	
	data_class = {
	    'bp': 'Biological Process', 
	    'cc': 'Cellular Component', 
	    'mf': 'Molecular Function'
	}
		
	# PDF prefix file name
	pdf_filename = '{0}_'.format(mirna_id)    
	
	if vtype == 'Substitution':
	    vpos = re.sub(r'[^0-9]+', '', request.form.getlist('vpos')[0]) if request.form.getlist('vpos')[0] else ''
	    vref = re.sub(r'[^AGCUT\-\+]+', '', request.form.getlist('vref')[0]) if request.form.getlist('vref')[0] else ''
	    vmod = re.sub(r'[^AGCUT\-\+]+', '', request.form.getlist('vmod')[0]) if request.form.getlist('vmod')[0] else ''
	    
	    if int(vpos) == -1 or vref == '' or vmod == '':
		return jsonify({ 'error_msg': 'No data available' }), 405
	    for s in stats['variations']['substitution']:
		if s['position'] == int(vpos) and \
		   s['ref'] == vref and \
		   s['mod'] == vmod:
		    section_title = '{0} (position: {1}, ref: {2}, mod: {3})'.format(
		        section_title, vpos, vref, vmod)
		    target = 'nt_var'
		    pdf_filename = '{0}_nt_var_{1}_{2}_{3}_{4}_'.format(
			pdf_filename,
			vpos,
			vref,
			vmod,
		        resource
		    )		
		    # For each GO class: bp, cc, mf
		    for c, c_name in data_class.iteritems():
			if c in s:
			    label = '{0}_his'.format(c)
			    x_labels = []
			    wt_data = []
			    var_data = []
			    for go in s[label]:
				if go['go_term'] in data_details:
				    lab = data_details[go['go_term']]['name']
				    x_labels.append(lab)
				    wt_data.append(float(go['wt']))
				    var_data.append(float(go['var']))
			    hist_data.append({
				'namespace': c_name,
				'x_labels': x_labels,
				'wt_data': wt_data,
				'var_data': var_data,
				'label': c,
				'threshold': [ stats['threshold'] for j in xrange(0, len(x_labels)) ]
			    })
		    break
		    
	elif vtype == 'isomiR':
	    v3p = re.sub(r'[^\-\+0-9]+', '', request.form.getlist('v3p')[0]) if request.form.getlist('v3p')[0] else ''
	    v5p = re.sub(r'[^\-\+0-9]+', '', request.form.getlist('v5p')[0]) if request.form.getlist('v5p')[0] else ''
	    
	    if v3p == '' or v5p == '':
		return jsonify({ 'error_msg': 'No data available' }), 405
	    
	    for i in stats['variations']['isomir']:
		if i['5p'] == v5p and i['3p'] == v3p:
		    section_title = '{0} (5p: {1}, 3p: {2})'.format(
		        section_title, v5p, v3p)		    
		    target = 'isomir'
		    pdf_filename = '{0}_isomir_{1}_{2}_{3}_'.format(
			pdf_filename,
			v5p,
			v3p,
		        resource
		    )		
		    # For each GO class: bp, cc, mf
		    for c, c_name in data_class.iteritems():
			if c in i:
			    label = '{0}_his'.format(c)
			    x_labels = []
			    wt_data = []
			    var_data = []
			    for go in i[label]:
				# Chart X label
				if go['go_term'] in data_details:
				    lab = data_details[go['go_term']]['name']
				    x_labels.append(lab)
				    wt_data.append(float(go['wt']))
				    var_data.append(float(go['var']))
			    hist_data.append({
				'namespace': c_name,
				'x_labels': x_labels,
				'wt_data': wt_data,
				'var_data': var_data,
				'label': c,
				'threshold': [ stats['threshold'] for j in xrange(0, len(x_labels)) ]
			    })
		    break
		
	elif vtype == 'Wild Type':
	    target = 'wildtype'
	    pdf_filename = '{0}_wildtype_{1}_'.format(pdf_filename, resource)		
	    # For each GO class: bp, cc, mf
	    for c, c_name in data_class.iteritems():
		if c in stats:
		    label = '{0}'.format(c)
		    x_labels = []
		    wt_data_dict = {}
		    wt_data = []
		    var_data = []
		    for go_term,details in stats[label].iteritems():
			# Chart X label
			if go_term in data_details:
			    lab = data_details[go_term]['name']			    
			    wt_data_dict[lab] = float(details['-Log(pvalue)'])
		    
		    for key, value in sorted(wt_data_dict.iteritems(), key=lambda (k,v): (v,k), reverse=True):
			x_labels.append(key)
			wt_data.append(value)
			
		    hist_data.append({
	                'namespace': c_name,
	                'x_labels': x_labels,
	                'wt_data': wt_data,
		        'var_data': var_data,
	                'label': c,
	                'threshold': [ stats['threshold'] for j in xrange(0, len(x_labels)) ]
	            })
    
    elif analysis_type_f == 'comparison':
	
	if 'list_id' not in data or \
	   'resource' not in data:
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	list_id = re.sub(r'[^0-9a-zA-Z\-]+', '', request.form.getlist('list_id')[0]) if request.form.getlist('list_id')[0] else ''
	resource = re.sub(r'[^a-zA-Z]+', '', request.form.getlist('resource')[0]) if request.form.getlist('resource')[0] else ''
	
	if list_id == '' or \
	   list_id not in ['list1', 'list2']:
	    return jsonify({ 'error_msg': 'No data available' }), 405    
	
	# GO term details
	data_details = {}
	if not os.path.isfile(res.GO_PATH):
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	with open(res.GO_PATH, 'r') as f:
	    data_details = json.load(f)
		
	enrich_data = {}
	if not os.path.isfile(res.ENRICH_COMP_JSON_PATH):
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	with open(res.ENRICH_COMP_JSON_PATH, 'r') as f:
	    enrich_data = json.load(f)	
	    
	if not enrich_data:
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	if list_id not in enrich_data[resource]:
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	stats = enrich_data[resource][list_id]
	
	data_class = {
	    'bp': 'Biological Process', 
	    'cc': 'Cellular Component', 
	    'mf': 'Molecular Function'
	}
		
	# PDF prefix file name
	pdf_filename = '{0}_{1}_{2}_'.format(enrich_data[resource]['list1']['entity'], 
	                                     stats['entity'],
	                                     resource)	
	
	entity_label1 = enrich_data[resource]['list1']['entity']
	entity_label2 = stats['entity']
	
	# For each GO class: bp, cc, mf
	for c, c_name in data_class.iteritems():
	    if c in stats:
		label = '{0}_his'.format(c)
		x_labels = []
		wt_data = []
		var_data = []
		for go in stats[label]:
		    # Chart X label
		    if go['go_term'] in data_details:
			lab = data_details[go['go_term']]['name']
			x_labels.append(lab)
			wt_data.append(float(go['wt']))
			var_data.append(float(go['var']))
		hist_data.append({
                    'namespace': c_name,
                    'x_labels': x_labels,
                    'wt_data': wt_data,
                    'var_data': var_data,
                    'label': c,
                    'threshold': [ stats['threshold'] for j in xrange(0, len(x_labels)) ]
                })	
	
    """Response"""    
    try:
        return render_template('enrichment_analysis/go_histogram.html',
	                       hist_data = hist_data,
	                       target = target,
	                       comparison = False,
	                       pdf_filename = pdf_filename,
	                       entity_label1 = entity_label1,
	                       entity_label2 = entity_label2,
	                       section_title = section_title)
    except TemplateNotFound:
        return jsonify({ 'error_msg': 'Not Found' }), 404
    
    
@enrichment_analysis.route('/enrichment_analysis_variation_comp_histograms', methods=['POST'])
def enrichment_analysis_variation_comp_histograms():
    
    data = request.form.to_dict()
    analysis_type = request.form.getlist('analysis_type')[0]
    
    # ==========================================================================
    # Input filtering 
    # ==========================================================================
    analysis_type_f = re.sub(r'[^a-z]+', '', analysis_type)
    
    if analysis_type_f not in ['comparison', 'single']:
	return jsonify({ 'error_msg': 'No data available' }), 405
    
    hist_data = []
    target = ''
    pdf_filename = ''
    entity_label1 = ''
    entity_label2 = '' 
    section_title = ''
    
    if analysis_type_f == 'single':
	    
	entity_label1 = 'Wild Type'
	entity_label2 = 'Variation'	
    
	if 'mirna_id' not in data or \
           'vtype' not in data or \
	   'resource' not in data:
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	mirna_id = re.sub(r'[^0-9a-zA-Z\-]+', '', request.form.getlist('mirna_id')[0]) if request.form.getlist('mirna_id')[0] else '' 
	vtype = re.sub(r'[^\sa-zA-Z]+', '', request.form.getlist('vtype')[0]) if request.form.getlist('vtype')[0] else ''
	resource = re.sub(r'[^a-zA-Z]+', '', request.form.getlist('resource')[0]) if request.form.getlist('resource')[0] else ''
	
	v3p = re.sub(r'[^\-\+0-9]+', '', request.form.getlist('v3p')[0]) if request.form.getlist('v3p')[0] else ''
	v5p = re.sub(r'[^\-\+0-9]+', '', request.form.getlist('v5p')[0]) if request.form.getlist('v5p')[0] else ''
	
	vpos = re.sub(r'[^0-9]+', '', request.form.getlist('vpos')[0]) if request.form.getlist('vpos')[0] else ''
	vref = re.sub(r'[^AGCUT\-\+]+', '', request.form.getlist('vref')[0]) if request.form.getlist('vref')[0] else ''
	vmod = re.sub(r'[^AGCUT\-\+]+', '', request.form.getlist('vmod')[0]) if request.form.getlist('vmod')[0] else ''
	
	if mirna_id == '' or \
           vtype == '' or \
           vtype not in ['isomiR', 'Substitution'] or \
	   resource not in ['go']:
	    return jsonify({ 'error_msg': 'No data available' }), 405    
	
	section_title = '{0}: {1}'.format(vtype, mirna_id)
	enrich_data = {}
	if not os.path.isfile(res.ENRICH_JSON_PATH):
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	with open(res.ENRICH_JSON_PATH, 'r') as f:
	    enrich_data = json.load(f)	
 
	if not enrich_data:
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	# GO term details
	data_details = {}
	if not os.path.isfile(res.GO_PATH):
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	with open(res.GO_PATH, 'r') as f:
	    data_details = json.load(f)
	
	if mirna_id not in enrich_data[resource]:
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	stats = enrich_data[resource][mirna_id]
	
	data_class = {
            'bp': 'Biological Process', 
            'cc': 'Cellular Component', 
            'mf': 'Molecular Function'
        }	
    
	# PDF prefix file name
	pdf_filename = '{0}_'.format(mirna_id)
	
	if vtype == 'Substitution':
	    if int(vpos) == -1 or vref == '' or vmod == '':
		return jsonify({ 'error_msg': 'No data available' }), 405
	    for s in stats['variations']['substitution']:
		if s['position'] == int(vpos) and \
		   s['ref'] == vref and \
		   s['mod'] == vmod:
		    section_title = '{0} (position: {1}, ref: {2}, mod: {3})'.format(
		        section_title, vpos, vref, vmod)		    
		    target = 'nt_var'
		    pdf_filename = '{0}_nt_var_{1}_{2}_{3}_{4}_'.format(
			pdf_filename,
			vpos,
			vref,
			vmod,
		        resource
		    )
		    # For each GO class: bp, cc, mf
		    for c, c_name in data_class.iteritems():
			if c in s:
			    label = 'comp_{0}_his'.format(c)
			    x_labels = []
			    wt_data = []
			    var_data = []
			    inter_data = []
			    for go in s[label]:
				# Chart X label
				if go['go_term'] in data_details:
				    lab = data_details[go['go_term']]['name']
				    x_labels.append(lab)
				    wt_data.append(float(go['wt_only']))
				    var_data.append(float(go['variation_only']))
				    inter_data.append(float(go['intersection']))
			    hist_data.append({
				'namespace': c_name,
				'x_labels': x_labels,
				'wt_data': wt_data,
				'var_data': var_data,
				'inter_data': inter_data,
				'label': c,
				'threshold': [ stats['threshold'] for j in xrange(0, len(x_labels)) ]
			    })
		    break
		    
	elif vtype == 'isomiR':
	    if v3p == '' or v5p == '':
		return jsonify({ 'error_msg': 'No data available' }), 405
	    for i in stats['variations']['isomir']:
		if i['5p'] == v5p and i['3p'] == v3p:
		    section_title = '{0} (5p: {1}, 3p: {2})'.format(
		        section_title, v5p, v3p)		    
		    target = 'isomir'
		    pdf_filename = '{0}_isomir_{1}_{2}_{3}_'.format(
			pdf_filename,
			v5p,
			v3p,
		        resource
		    )		
		    # For each GO class: bp, cc, mf
		    for c, c_name in data_class.iteritems():
			if c in i:
			    label = 'comp_{0}_his'.format(c)
			    x_labels = []
			    wt_data = []
			    var_data = []
			    inter_data = []
			    for go in i[label]:
				# Chart X label
				if go['go_term'] in data_details:
				    lab = data_details[go['go_term']]['name']
				    x_labels.append(lab)
				    wt_data.append(float(go['wt_only']))
				    var_data.append(float(go['variation_only']))
				    inter_data.append(float(go['intersection']))
			    hist_data.append({
				'namespace': c_name,
				'x_labels': x_labels,
				'wt_data': wt_data,
				'var_data': var_data,
				'inter_data': inter_data,
				'label': c,
				'threshold': [ stats['threshold'] for j in xrange(0, len(x_labels)) ]
			    })
		    break
		
    elif analysis_type_f == 'comparison':

	if 'list_id' not in data or \
	   'resource' not in data:
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	list_id = re.sub(r'[^0-9a-zA-Z\-]+', '', request.form.getlist('list_id')[0]) if request.form.getlist('list_id')[0] else ''
	resource = re.sub(r'[^a-zA-Z]+', '', request.form.getlist('resource')[0]) if request.form.getlist('resource')[0] else ''
	
	if list_id == '' or \
	   list_id not in ['list1', 'list2'] or \
	   resource not in ['go']:
	    return jsonify({ 'error_msg': 'No data available' }), 405    
	
	# GO term details
	data_details = {}
	if not os.path.isfile(res.GO_PATH):
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	with open(res.GO_PATH, 'r') as f:
	    data_details = json.load(f)
		
	enrich_data = {}
	if not os.path.isfile(res.ENRICH_COMP_JSON_PATH):
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	with open(res.ENRICH_COMP_JSON_PATH, 'r') as f:
	    enrich_data = json.load(f)	
	    
	if not enrich_data:
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	if list_id not in enrich_data[resource]:
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	stats = enrich_data[resource][list_id]
	
	go_class = {
	    'bp': 'Biological Process', 
	    'cc': 'Cellular Component', 
	    'mf': 'Molecular Function'
	}
		
	# PDF prefix file name
	pdf_filename = '{0}_{1}_{2}_'.format(enrich_data[resource]['list1']['entity'], 
	                                     stats['entity'],
	                                     resource)	
	
	entity_label1 = enrich_data[resource]['list1']['entity']
	entity_label2 = stats['entity']
	
	# For each GO class: bp, cc, mf
	for c, c_name in go_class.iteritems():
	    if c in stats:
		label = 'comp_{0}_his'.format(c)
		x_labels = []
		wt_data = []
		var_data = []
		inter_data = []
		for go in stats[label]:
		    # Chart X label
		    if go['go_term'] in data_details:
			lab = data_details[go['go_term']]['name']
			x_labels.append(lab)
			wt_data.append(float(go['wt_only']))
			var_data.append(float(go['variation_only']))
			inter_data.append(float(go['intersection']))
		hist_data.append({
                    'namespace': c_name,
		    'x_labels': x_labels,
		    'wt_data': wt_data,
		    'var_data': var_data,
		    'inter_data': inter_data,
		    'label': c,
		    'threshold': [ stats['threshold'] for j in xrange(0, len(x_labels)) ]
                })
    
    """Response"""    
    try:
        return render_template('enrichment_analysis/go_histogram.html',
	                       hist_data = hist_data,
	                       target = target,
	                       comparison = True,
	                       pdf_filename = pdf_filename,
	                       entity_label1 = entity_label1,
	                       entity_label2 = entity_label2,
	                       section_title = section_title)
    except TemplateNotFound:
        return jsonify({ 'error_msg': 'Not Found' }), 404
    
    
@enrichment_analysis.route('/enrichment_analysis_kill_all', methods=['POST'])
def enrichment_analysis_kill_all():
    
    if request.method != 'POST':
        return jsonify({ 'error_msg': 'No data available' }), 405
    
    analysis_type = request.form.getlist('analysis_type')[0]
    
    # ==========================================================================
    # Input filtering 
    # ==========================================================================
    analysis_type_f = re.sub(r'[^a-z]+', '', analysis_type)
    
    if analysis_type_f not in ['comparison', 'single']:
	return jsonify({ 'error_msg': 'No data available' }), 405
    
    f_path = res.JOB_ENRICH_PATH
    if analysis_type_f == 'comparison':  
	f_path = res.JOB_ENRICH_COMP_PATH	
	
    job_data = {}
    if os.path.isfile(f_path):
	with open(f_path, 'r') as f:
	    job_data = json.load(f)	
   
    if 'job_pid' in job_data:
	pid = int(job_data['job_pid'])  
	if pid:
	    # If process is not a process group
	    os.system('pkill -TERM -P {pid}'.format(pid = pid))
	
	    parent = psutil.Process(pid)
	    for child in parent.children(recursive=True):  # or parent.children() for recursive=False
		child.kill()
	    parent.kill()
	
    # Remove zombie processes
    process = subprocess.Popen(['bash','-c', '/opt/kill_zombies.sh'], stdout=subprocess.PIPE)
    output, error = process.communicate()
    
    remove_file(f_path)
        
    return jsonify({ 'status': 'done' }), 200


'''
================================================================================
================================================================================
================================================================================
Individual Functional Enrichment Analysis
================================================================================
================================================================================
================================================================================
'''
@enrichment_analysis.route('/exec_functional_enrichment_comp_analysis', methods=['POST'])
def exec_functional_enrichment_comp_analysis():
    
    if request.method != 'POST':
        return jsonify({ 'error_msg': 'No data available' }), 405  
    
    if os.path.isfile(res.JOB_PATH) or \
       os.path.isfile(res.JOB_ENRICH_PATH) or \
       os.path.isfile(res.JOB_ENRICH_COMP_PATH):
	return jsonify({ 'error_msg': 'There is already a running job.' }), 405    
    
    cores = request.form.getlist('cores')[0]
    targets_l1 = request.form.getlist('l1')[0]
    targets_l2 = request.form.getlist('l2')[0]
    name_l1 = request.form.getlist('name_l1')[0]
    name_l2 = request.form.getlist('name_l2')[0] 
    pvalue = request.form.getlist('pvalue')[0]
    pvalue_type = request.form.getlist('pvalue_type')[0]
    enrich_opts = request.form.getlist('enrich_opts[]')
    
    # ==========================================================================
    # Input filtering 
    # ==========================================================================
    cores_f = int(re.sub(r'[^0-9]+', '', cores))
    targets_l1_f = re.sub(r'[^0-9a-zA-Z_\s\.]+', '', targets_l1)
    targets_l2_f = re.sub(r'[^0-9a-zA-Z_\s\.]+', '', targets_l2)
    name_l1_f = re.sub(r'\s+', '_', re.sub(r'[^0-9a-zA-Z_\s\-]+', '', name_l1))
    name_l2_f = re.sub(r'\s+', '_', re.sub(r'[^0-9a-zA-Z_\s\-]+', '', name_l2))  
    pvalue_f = float(re.sub(r'[^0-9\.]+', '', pvalue))
    pvalue_type_f = re.sub(r'[^a-z]+', '', pvalue_type)  
    enrich_opts_f = []
    
    for el in enrich_opts:
        f_el = re.sub(r'[^a-z]', '', el)
        if f_el in ['go', 'kegg', 'reactome']:
            enrich_opts_f.append(f_el) 
        	
    if not isinstance(cores_f, int) or \
       not isinstance(pvalue_f, float) or \
       pvalue_type_f not in ['raw', 'adj'] or \
       targets_l1_f == '' or \
       name_l1_f == '' or \
       len(enrich_opts_f) == 0:
        return jsonify({ 'error_msg': 'No data available' }), 405
    
    if name_l1_f == name_l2_f:
        return jsonify({ 'error_msg': 'Both lists have the same name' }), 405
    
    if targets_l1_f == targets_l2_f:
        return jsonify({ 'error_msg': 'Both lists have the targets set' }), 405 
    
    if name_l2_f != '' and targets_l2_f == '':
        return jsonify({ 'error_msg': 'No targets list has been specified for the second list' }), 405 
    
    if name_l2_f == '' and targets_l2_f != '':
        return jsonify({ 'error_msg': 'No targets list has been specified for the second list' }), 405    
    
    command = "/usr/bin/python /app/iso/run_comp.py --cores {0} --list1 {1} --list2 {2} --name1 {3} --name2 {4} --pvalue {5} --pvalue_type {6} --enrichtype {7}".format(
        cores_f, ','.join(targets_l1_f.split()), ','.join(targets_l2_f.split()), 
        name_l1_f, name_l2_f, pvalue_f, pvalue_type_f, ','.join(enrich_opts_f))
    
    if targets_l2_f == '' or name_l2_f == '':
	command = "/usr/bin/python /app/iso/run_comp.py --cores {0} --list1 {1} --name1 {2} --pvalue {3} --pvalue_type {4} --enrichtype {5}".format(
	    cores_f, ','.join(targets_l1_f.split()), name_l1_f, pvalue_f, 
	    pvalue_type_f, ','.join(enrich_opts_f))
	
    try:
	task = celery_run_command.delay(command.split()) 
	task_id = task.id
	task_status = celery_run_command.AsyncResult(task_id)
	job_details = task.get()
	
	job = {
	    'task_id': task_id, 
	    'job_pid': job_details['pid'],
	    'job_status': job_details['status'],
	    'job_msg': job_details['msg'],
	    'command': command
	}
	
	with open(res.JOB_ENRICH_COMP_PATH, 'w') as f:
	    f.write(json.dumps(job, indent=4))	
	
	return jsonify(job), 202

    except ValidationError as e:
	return jsonify({ 'error_msg': e.messages }), 405
    
    
@enrichment_analysis.route('/reload_enrichment_analysis_results', methods=['POST'])
def reload_enrichment_analysis_results():

    if request.method != 'POST':
        return jsonify({ 'error_msg': 'No data available' }), 405
	   
    """Response"""    
    try:
        return render_template('enrichment_analysis/enrich_reload_results.html')
    except TemplateNotFound:
        return jsonify({ 'error_msg': 'Not Found' }), 404
    
    
@enrichment_analysis.route('/pathway_enrichment_analysis_variation_statistics', methods=['POST'])   
def pathway_enrichment_analysis_variation_statistics():
    
    if request.method != 'POST':
        return jsonify({ 'error_msg': 'No data available' }), 405
    
    data = request.form.to_dict()
    analysis_type = request.form.getlist('analysis_type')[0]
    
    # ==========================================================================
    # Input filtering 
    # ==========================================================================
    analysis_type_f = re.sub(r'[^a-z]+', '', analysis_type)
    
    if analysis_type_f not in ['comparison', 'single']:
	return jsonify({ 'error_msg': 'No data available' }), 405
    
    table_data = []
    filename = ''
    
    if analysis_type_f == 'single':	
	
	if 'mirna_id' not in data or \
	   'vtype' not in data or \
	   'data_tab' not in data or \
	   'resource' not in data:
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	mirna_id = re.sub(r'[^0-9a-zA-Z\-]+', '', request.form.getlist('mirna_id')[0]) if request.form.getlist('mirna_id')[0] else '' 
	vtype = re.sub(r'[^\sa-zA-Z]+', '', request.form.getlist('vtype')[0]) if request.form.getlist('vtype')[0] else ''
	data_tab = re.sub(r'[^\sa-zA-Z\_]+', '', request.form.getlist('data_tab')[0]) if request.form.getlist('data_tab')[0] else ''
	resource = re.sub(r'[^a-zA-Z]+', '', request.form.getlist('resource')[0]) if request.form.getlist('resource')[0] else ''
	data_tab = data_tab.replace('enrich_row_', '')
	
	if mirna_id == '' or \
	   vtype == '' or \
	   vtype not in ['Wild Type', 'isomiR', 'Substitution'] or \
	   data_tab not in [
	       'pathway', 'pathway_wt_o', 'pathway_ins', 'pathway_var_o'
	    ] or resource not in ['kegg', 'reactome']:
	    return jsonify({ 'error_msg': 'No data available' }), 405    
	
	enrich_data = {}
	if not os.path.isfile(res.ENRICH_JSON_PATH):
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	with open(res.ENRICH_JSON_PATH, 'r') as f:
	    enrich_data = json.load(f)	
	    
	if not enrich_data:
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	if mirna_id not in enrich_data[resource]:
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	stats = enrich_data[resource][mirna_id]
	
	# Check the selected data category
	data_cat = ''    
	data_label = data_tab.split('_')
	
	if len(data_label) == 1:
	    pass
	elif len(data_label) == 2:
	    data_cat = 'intersection'
	elif len(data_label) == 3:
	    data_cat = 'wt_only'  
	    if data_label[1] == 'var':
		data_cat = 'variation_only'
	
	if vtype == 'Substitution':
	    vpos = re.sub(r'[^0-9]+', '', request.form.getlist('vpos')[0]) if request.form.getlist('vpos')[0] else ''
	    vref = re.sub(r'[^AGCUT\-\+]+', '', request.form.getlist('vref')[0]) if request.form.getlist('vref')[0] else ''
	    vmod = re.sub(r'[^AGCUT\-\+]+', '', request.form.getlist('vmod')[0]) if request.form.getlist('vmod')[0] else ''
	    
	    if int(vpos) == -1 or vref == '' or vmod == '':
		return jsonify({ 'error_msg': 'No data available' }), 405
	    for s in stats['variations']['substitution']:
		if s['position'] == int(vpos) and \
	           s['ref'] == vref and \
	           s['mod'] == vmod:			
		    if data_cat != '':
			if data_cat in s:
			    table_data = s[data_cat]
			    filename = '{0}_nt_var_{1}_{2}_{3}_{4}_{5}'.format(mirna_id, 
			                                                       s['position'],
			                                                       s['ref'],
			                                                       s['mod'],
			                                                       data_cat,
			                                                       resource)				
		    else:
			table_data = s['pathways']
			filename = '{0}_nt_var_{1}_{2}_{3}_{4}'.format(mirna_id, 
			                                               s['position'],
			                                               s['ref'],
			                                               s['mod'],
			                                               resource)			    
		    break
			
	elif vtype == 'isomiR':
	    v3p = re.sub(r'[^\-\+0-9]+', '', request.form.getlist('v3p')[0]) if request.form.getlist('v3p')[0] else ''
	    v5p = re.sub(r'[^\-\+0-9]+', '', request.form.getlist('v5p')[0]) if request.form.getlist('v5p')[0] else ''
	    
	    if v3p == '' or v5p == '':
		return jsonify({ 'error_msg': 'No data available' }), 405
	    for i in stats['variations']['isomir']:
		if i['5p'] == v5p and i['3p'] == v3p:			
		    if data_cat != '':
			if data_cat in i:
			    table_data = i[data_cat]
			    filename = '{0}_isomir_{1}_{2}_{3}_{4}'.format(mirna_id, 
			                                                   i['5p'],
			                                                   i['3p'],
			                                                   data_cat,
			                                                   resource)				
		    else:
			table_data = i['pathways']
			filename = '{0}_isomir_{1}_{2}_{3}'.format(mirna_id, 
			                                           i['5p'],
			                                           i['3p'],
			                                           resource)			    
		    break
			
	elif vtype == 'Wild Type':
		filename = '{0}_wild_type_{1}'.format(mirna_id, resource)
		table_data = stats['pathways']
    
    elif analysis_type_f == 'comparison':

	if 'list_id' not in data or \
	   'data_tab' not in data or \
	   'resource' not in data:
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	list_id = re.sub(r'[^0-9a-zA-Z\-]+', '', request.form.getlist('list_id')[0]) if request.form.getlist('list_id')[0] else '' 
	resource = re.sub(r'[^a-zA-Z]+', '', request.form.getlist('resource')[0]) if request.form.getlist('resource')[0] else ''
	data_tab = re.sub(r'[^\sa-zA-Z\_]+', '', request.form.getlist('data_tab')[0]) if request.form.getlist('data_tab')[0] else ''
	data_tab = data_tab.replace('enrich_row_', '')
	
	if list_id == '' or \
	   list_id not in ['list1', 'list2'] or \
	   data_tab not in [
	       'pathway', 'pathway_wt_o', 'pathway_ins', 'pathway_var_o'
	    ] or resource not in ['kegg', 'reactome']:
	    return jsonify({ 'error_msg': 'No data available' }), 405    
	
	enrich_data = {}
	if not os.path.isfile(res.ENRICH_COMP_JSON_PATH):
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	with open(res.ENRICH_COMP_JSON_PATH, 'r') as f:
	    enrich_data = json.load(f)	
	    
	if not enrich_data:
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	if list_id not in enrich_data[resource]:
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	stats = enrich_data[resource][list_id]	
	
	# Check the selected GO term category and class
	data_cat = ''
	data_label = data_tab.split('_')
	
	if len(data_label) == 1:
	    pass
	elif len(data_label) == 2:
	    data_cat = 'intersection'
	    data_cat_f = 'intersection'
	elif len(data_label) == 3:
	    data_cat = 'wt_only'  
	    data_cat_f = 'list1_only'
	    if data_label[1] == 'var':
		data_cat = 'variation_only'
		data_cat_f = 'list2_only'
		
	if data_cat != '':
	    if data_cat in stats:
		table_data = stats[data_cat]
		filename = '{0}_{1}_{2}'.format(stats['entity'], 
		                                data_cat_f,
		                                resource)
	else:
	    table_data = stats['pathways']	
	    filename = '{0}_{1}'.format(stats['entity'], resource)		
    
    """Response"""    
    try:
        return render_template('enrichment_analysis/pathway_statistics.html',
	                       table_data = table_data,
	                       filename = filename,
	                       resource = resource)
    except TemplateNotFound:
        return jsonify({ 'error_msg': 'Not Found' }), 404
    
    
    
@enrichment_analysis.route('/pathway_enrichment_analysis_variation_histograms', methods=['POST'])
def pathway_enrichment_analysis_variation_histograms():
    
    data = request.form.to_dict()
    analysis_type = request.form.getlist('analysis_type')[0]
    
    # ==========================================================================
    # Input filtering 
    # ==========================================================================
    analysis_type_f = re.sub(r'[^a-z]+', '', analysis_type)
    
    if analysis_type_f not in ['comparison', 'single']:
	return jsonify({ 'error_msg': 'No data available' }), 405
    
    hist_data = []
    target = ''
    pdf_filename = ''
    entity_label1 = ''
    entity_label2 = ''    
    
    section_title = ''
    
    if analysis_type_f == 'single':
	
	entity_label1 = 'Wild Type'
	entity_label2 = 'Variation'	
    
	if 'mirna_id' not in data or \
	   'vtype' not in data or \
	   'resource' not in data:
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	mirna_id = re.sub(r'[^0-9a-zA-Z\-]+', '', request.form.getlist('mirna_id')[0]) if request.form.getlist('mirna_id')[0] else '' 
	vtype = re.sub(r'[^\sa-zA-Z]+', '', request.form.getlist('vtype')[0]) if request.form.getlist('vtype')[0] else ''
	resource = re.sub(r'[^a-zA-Z]+', '', request.form.getlist('resource')[0]) if request.form.getlist('resource')[0] else ''
	
	if mirna_id == '' or \
	   vtype == '' or \
	   vtype not in ['Wild Type', 'isomiR', 'Substitution'] or \
	   resource not in ['kegg', 'reactome']:
	    return jsonify({ 'error_msg': 'No data available' }), 405    
	
	section_title = '{0}: {1}'.format(vtype, mirna_id)
	
	enrich_data = {}
	if not os.path.isfile(res.ENRICH_JSON_PATH):
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	with open(res.ENRICH_JSON_PATH, 'r') as f:
	    enrich_data = json.load(f)	
 
	if not enrich_data:
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	# Pathway details
	data_details = {}
	if not os.path.isfile(res.PATHWAY_PATH):
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	with open(res.PATHWAY_PATH, 'r') as f:
	    d = json.load(f)
	    if resource not in d:
		return jsonify({ 'error_msg': 'No data available' }), 405
	    data_details = d[resource]
	
	if mirna_id not in enrich_data[resource]:
	    return jsonify({ 'error_msg': 'No data available' }), 405

	stats = enrich_data[resource][mirna_id]
		
	# PDF prefix file name
	pdf_filename = '{0}_'.format(mirna_id)    
	
	if vtype == 'Substitution':
	    vpos = re.sub(r'[^0-9]+', '', request.form.getlist('vpos')[0]) if request.form.getlist('vpos')[0] else ''
	    vref = re.sub(r'[^AGCUT\-\+]+', '', request.form.getlist('vref')[0]) if request.form.getlist('vref')[0] else ''
	    vmod = re.sub(r'[^AGCUT\-\+]+', '', request.form.getlist('vmod')[0]) if request.form.getlist('vmod')[0] else ''
	    
	    if int(vpos) == -1 or vref == '' or vmod == '':
		return jsonify({ 'error_msg': 'No data available' }), 405
	    for s in stats['variations']['substitution']:
		if s['position'] == int(vpos) and \
		   s['ref'] == vref and \
		   s['mod'] == vmod:
		    section_title = '{0} (position: {1}, ref: {2}, mod: {3})'.format(
		        section_title, vpos, vref, vmod)		    
		    target = 'nt_var'
		    pdf_filename = '{0}_nt_var_{1}_{2}_{3}_'.format(
			pdf_filename,
			vpos,
			vref,
			vmod
		    )		
		    
		    x_labels = []
		    wt_data = []
		    var_data = []
		    for pathway in s['his']:
			if pathway['pathway_id'] in data_details:
			    lab = data_details[pathway['pathway_id']]['pathway_name']
			    x_labels.append(lab)
			    wt_data.append(float(pathway['wt']))
			    var_data.append(float(pathway['var']))
		    hist_data.append({
		        'label': resource,
	                'x_labels': x_labels,
	                'wt_data': wt_data,
	                'var_data': var_data,
	                'threshold': [ stats['threshold'] for j in xrange(0, len(x_labels)) ]
	            })
		    break
		    
	elif vtype == 'isomiR':
	    v3p = re.sub(r'[^\-\+0-9]+', '', request.form.getlist('v3p')[0]) if request.form.getlist('v3p')[0] else ''
	    v5p = re.sub(r'[^\-\+0-9]+', '', request.form.getlist('v5p')[0]) if request.form.getlist('v5p')[0] else ''
	    
	    if v3p == '' or v5p == '':
		return jsonify({ 'error_msg': 'No data available' }), 405
	    
	    for i in stats['variations']['isomir']:
		if i['5p'] == v5p and i['3p'] == v3p:
		    section_title = '{0} (5p: {1}, 3p: {2})'.format(
		        section_title, v5p, v3p)		    
		    target = 'isomir'
		    pdf_filename = '{0}_isomir_{1}_{2}_'.format(
			pdf_filename,
			v5p,
			v3p
		    )	
		    
		    x_labels = []
		    wt_data = []
		    var_data = []
		    for pathway in i['his']:
			if pathway['pathway_id'] in data_details:
			    lab = data_details[pathway['pathway_id']]['pathway_name']
			    x_labels.append(lab)
			    wt_data.append(float(pathway['wt']))
			    var_data.append(float(pathway['var']))
		    hist_data.append({
		        'label': resource,
	                'x_labels': x_labels,
	                'wt_data': wt_data,
	                'var_data': var_data,
	                'threshold': [ stats['threshold'] for j in xrange(0, len(x_labels)) ]
	            })
		    break
		
	elif vtype == 'Wild Type':
	    target = 'wildtype'
	    pdf_filename = '{0}_wildtype_'.format(pdf_filename)		
	   
	    label = 'pathways'
	    x_labels = []
	    wt_data_dict = {}
	    wt_data = []
	    var_data = []
	    for pathway_id,details in stats[label].iteritems():
		# Chart X label
		if pathway_id in data_details:
		    lab = data_details[pathway_id]['pathway_name']
		    wt_data_dict[lab] = float(details['-Log(pvalue)'])
		    
	    for key, value in sorted(wt_data_dict.iteritems(), key=lambda (k,v): (v,k), reverse=True):
		x_labels.append(key)
		wt_data.append(value)
	
	    hist_data.append({
                'x_labels': x_labels,
                'wt_data': wt_data,
                'var_data': var_data,
                'label': resource,
                'threshold': [ stats['threshold'] for j in xrange(0, len(x_labels)) ]
            })
    
    elif analysis_type_f == 'comparison':
	if 'list_id' not in data or \
	   'resource' not in data:
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	resource = re.sub(r'[^a-zA-Z]+', '', request.form.getlist('resource')[0]) if request.form.getlist('resource')[0] else ''	
	list_id = re.sub(r'[^0-9a-zA-Z\-]+', '', request.form.getlist('list_id')[0]) if request.form.getlist('list_id')[0] else '' 
	
	if list_id == '' or \
	   list_id not in ['list1', 'list2'] or \
	   resource not in ['kegg', 'reactome']:
	    return jsonify({ 'error_msg': 'No data available' }), 405    
	
	enrich_data = {}
	if not os.path.isfile(res.ENRICH_COMP_JSON_PATH):
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	with open(res.ENRICH_COMP_JSON_PATH, 'r') as f:
	    enrich_data = json.load(f)	
 
	if not enrich_data:
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	# Pathway details
	data_details = {}
	if not os.path.isfile(res.PATHWAY_PATH):
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	with open(res.PATHWAY_PATH, 'r') as f:
	    d = json.load(f)
	    if resource not in d:
		return jsonify({ 'error_msg': 'No data available' }), 405
	    data_details = d[resource]
	
	if list_id not in enrich_data[resource]:
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	stats = enrich_data[resource][list_id]
		
	# PDF prefix file name
	pdf_filename = '{0}_{1}_'.format(enrich_data[resource]['list1']['entity'], 
	                                    stats['entity'])	
	
	entity_label1 = enrich_data[resource]['list1']['entity']
	entity_label2 = stats['entity']
	
	x_labels = []
	wt_data = []
	var_data = []
	for pathway in stats['his']:
	    if pathway['pathway_id'] in data_details:
		lab = data_details[pathway['pathway_id']]['pathway_name']
		x_labels.append(lab)
		wt_data.append(float(pathway['wt']))
		var_data.append(float(pathway['var']))
	hist_data.append({
            'x_labels': x_labels,
            'wt_data': wt_data,
            'var_data': var_data,
            'label': resource,
            'threshold': [ stats['threshold'] for j in xrange(0, len(x_labels)) ]
        })	
	
    """Response"""    
    try:
        return render_template('enrichment_analysis/pathway_histogram.html',
	                       hist_data = hist_data,
	                       target = target,
	                       comparison = False,
	                       pdf_filename = pdf_filename,
	                       entity_label1 = entity_label1,
	                       entity_label2 = entity_label2,
	                       resource = resource,
	                       section_title = section_title)
    except TemplateNotFound:
        return jsonify({ 'error_msg': 'Not Found' }), 404
    
    
@enrichment_analysis.route('/pathway_enrichment_analysis_variation_comp_histograms', methods=['POST'])
def pathway_enrichment_analysis_variation_comp_histograms():
    
    data = request.form.to_dict()
    analysis_type = request.form.getlist('analysis_type')[0]
    
    # ==========================================================================
    # Input filtering 
    # ==========================================================================
    analysis_type_f = re.sub(r'[^a-z]+', '', analysis_type)
    
    if analysis_type_f not in ['comparison', 'single']:
	return jsonify({ 'error_msg': 'No data available' }), 405
    
    hist_data = []
    target = ''
    pdf_filename = ''
    entity_label1 = ''
    entity_label2 = ''  
    section_title = ''
    
    if analysis_type_f == 'single':
	    
	entity_label1 = 'Wild Type'
	entity_label2 = 'Variation'	
    
	if 'mirna_id' not in data or \
           'vtype' not in data:
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	mirna_id = re.sub(r'[^0-9a-zA-Z\-]+', '', request.form.getlist('mirna_id')[0]) if request.form.getlist('mirna_id')[0] else '' 
	vtype = re.sub(r'[^\sa-zA-Z]+', '', request.form.getlist('vtype')[0]) if request.form.getlist('vtype')[0] else ''
	resource = re.sub(r'[^a-zA-Z]+', '', request.form.getlist('resource')[0]) if request.form.getlist('resource')[0] else ''
	
	v3p = re.sub(r'[^\-\+0-9]+', '', request.form.getlist('v3p')[0]) if request.form.getlist('v3p')[0] else ''
	v5p = re.sub(r'[^\-\+0-9]+', '', request.form.getlist('v5p')[0]) if request.form.getlist('v5p')[0] else ''
	
	vpos = re.sub(r'[^0-9]+', '', request.form.getlist('vpos')[0]) if request.form.getlist('vpos')[0] else ''
	vref = re.sub(r'[^AGCUT\-\+]+', '', request.form.getlist('vref')[0]) if request.form.getlist('vref')[0] else ''
	vmod = re.sub(r'[^AGCUT\-\+]+', '', request.form.getlist('vmod')[0]) if request.form.getlist('vmod')[0] else ''
	
	if mirna_id == '' or \
	   vtype == '' or \
	   vtype not in ['isomiR', 'Substitution'] or \
	   resource not in ['kegg', 'reactome']:
	    return jsonify({ 'error_msg': 'No data available' }), 405    
	
	section_title = '{0}: {1}'.format(vtype, mirna_id)
	enrich_data = {}
	if not os.path.isfile(res.ENRICH_JSON_PATH):
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	with open(res.ENRICH_JSON_PATH, 'r') as f:
	    enrich_data = json.load(f)	
 
	if not enrich_data:
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	# Pathway details
	data_details = {}
	if not os.path.isfile(res.PATHWAY_PATH):
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	with open(res.PATHWAY_PATH, 'r') as f:
	    d = json.load(f)
	    if resource not in d:
		return jsonify({ 'error_msg': 'No data available' }), 405
	    data_details = d[resource]
	
	if mirna_id not in enrich_data[resource]:
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	stats = enrich_data[resource][mirna_id]
	
	# PDF prefix file name
	pdf_filename = '{0}_'.format(mirna_id)
	
	if vtype == 'Substitution':
	    if int(vpos) == -1 or vref == '' or vmod == '':
		return jsonify({ 'error_msg': 'No data available' }), 405
	    for s in stats['variations']['substitution']:
		if s['position'] == int(vpos) and \
		   s['ref'] == vref and \
		   s['mod'] == vmod:
		    section_title = '{0} (position: {1}, ref: {2}, mod: {3})'.format(
		        section_title, vpos, vref, vmod)		    
		    target = 'nt_var'
		    pdf_filename = '{0}_nt_var_{1}_{2}_{3}_'.format(
			pdf_filename,
			vpos,
			vref,
			vmod
		    )
			    
		    x_labels = []
		    wt_data = []
		    var_data = []
		    inter_data = []
		    for pathway in s['comp_his']:
			# Chart X label
			if pathway['pathway_id'] in data_details:
			    lab = data_details[pathway['pathway_id']]['pathway_name']
			    x_labels.append(lab)
			    wt_data.append(float(pathway['wt_only']))
			    var_data.append(float(pathway['variation_only']))
			    inter_data.append(float(pathway['intersection']))
		    hist_data.append({
	                'x_labels': x_labels,
	                'wt_data': wt_data,
	                'var_data': var_data,
	                'inter_data': inter_data,
	                'label': resource,
	                'threshold': [ stats['threshold'] for j in xrange(0, len(x_labels)) ]
	            })
		    break
		    
	elif vtype == 'isomiR':
	    if v3p == '' or v5p == '':
		return jsonify({ 'error_msg': 'No data available' }), 405
	    for i in stats['variations']['isomir']:
		if i['5p'] == v5p and i['3p'] == v3p:
		    section_title = '{0} (5p: {1}, 3p: {2})'.format(
		        section_title, v5p, v3p)		    
		    target = 'isomir'
		    pdf_filename = '{0}_isomir_{1}_{2}_'.format(
			pdf_filename,
			v5p,
			v3p
		    )		
		    
		    x_labels = []
		    wt_data = []
		    var_data = []
		    inter_data = []
		    for pathway in i['comp_his']:
			# Chart X label
			if pathway['pathway_id'] in data_details:
			    lab = data_details[pathway['pathway_id']]['pathway_name']
			    x_labels.append(lab)
			    wt_data.append(float(pathway['wt_only']))
			    var_data.append(float(pathway['variation_only']))
			    inter_data.append(float(pathway['intersection']))
		    hist_data.append({
	                'x_labels': x_labels,
	                'wt_data': wt_data,
	                'var_data': var_data,
	                'inter_data': inter_data,
	                'label': resource,
	                'threshold': [ stats['threshold'] for j in xrange(0, len(x_labels)) ]
	            })
		    break
		
    elif analysis_type_f == 'comparison':
	if 'list_id' not in data or \
	   'resource' not in data:
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	list_id = re.sub(r'[^0-9a-zA-Z\-]+', '', request.form.getlist('list_id')[0]) if request.form.getlist('list_id')[0] else '' 
	resource = re.sub(r'[^a-zA-Z]+', '', request.form.getlist('resource')[0]) if request.form.getlist('resource')[0] else ''
	
	if list_id == '' or \
	   list_id not in ['list1', 'list2'] or \
	   resource not in ['kegg', 'reactome']:
	    return jsonify({ 'error_msg': 'No data available' }), 405    
	
	enrich_data = {}
	if not os.path.isfile(res.ENRICH_COMP_JSON_PATH):
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	with open(res.ENRICH_COMP_JSON_PATH, 'r') as f:
	    enrich_data = json.load(f)	
 
	if not enrich_data:
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	# Pathway details
	data_details = {}
	if not os.path.isfile(res.PATHWAY_PATH):
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	with open(res.PATHWAY_PATH, 'r') as f:
	    d = json.load(f)
	    if resource not in d:
		return jsonify({ 'error_msg': 'No data available' }), 405
	    data_details = d[resource]
	
	if list_id not in enrich_data[resource]:
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	stats = enrich_data[resource][list_id]
		
	# PDF prefix file name
	pdf_filename = '{0}_{1}_'.format(enrich_data[resource]['list1']['entity'], 
	                                 stats['entity'])	
	
	entity_label1 = enrich_data[resource]['list1']['entity']
	entity_label2 = stats['entity']
	
	x_labels = []
	wt_data = []
	var_data = []
	inter_data = []
	for pathway in stats['comp_his']:
	    if pathway['pathway_id'] in data_details:
		lab = data_details[pathway['pathway_id']]['pathway_name']
		x_labels.append(lab)
		wt_data.append(float(pathway['wt_only']))
		var_data.append(float(pathway['variation_only']))
		inter_data.append(float(pathway['intersection']))
	hist_data.append({
            'x_labels': x_labels,
            'wt_data': wt_data,
            'var_data': var_data,
            'inter_data': inter_data,
            'label': resource,
            'threshold': [ stats['threshold'] for j in xrange(0, len(x_labels)) ]
        })
    
    """Response"""    
    try:
        return render_template('enrichment_analysis/pathway_histogram.html',
	                       hist_data = hist_data,
	                       target = target,
	                       comparison = True,
	                       pdf_filename = pdf_filename,
	                       entity_label1 = entity_label1,
	                       entity_label2 = entity_label2,
	                       section_title = section_title)
    except TemplateNotFound:
        return jsonify({ 'error_msg': 'Not Found' }), 404