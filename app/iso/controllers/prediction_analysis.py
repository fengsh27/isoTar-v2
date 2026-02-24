# -*- coding: UTF-8 -*-
# Import flask dependencies
from flask import (
    Blueprint, 
    render_template, 
    make_response,
    jsonify,
    request,
    Response,
    current_app
)
from jinja2 import TemplateNotFound
import requests
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

prediction_analysis = Blueprint('prediction_analysis', __name__)


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
	pass
    
    return retcode 

def remove_file(path):
    if os.path.isfile(path):
	try:
	    os.remove(path)
	except OSError as e:
	    if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
		pass
	
	
@prediction_analysis.route('/prediction_analysis')
def index_prediction_analysis():
    
    # Get the number of CPU cores
    #no_system_cores = psutil.cpu_count()
    no_system_cores = available_cpu_count()
    
    # Mature miRNAs (Mirbase v20)
    mat_mirnas_list = []
    if os.path.isfile(res.PRE_MAT_MIRNA_PATH):
        with open(res.PRE_MAT_MIRNA_PATH, 'r') as f:
            mirnas = json.load(f) 
            #mat_mirnas_list = list(map(str, mirnas.keys()))
	    for mirna_id,details in mirnas.iteritems():		
		for pre in details:		    
		    line = '%s [pre-miR: %s]' % (str(mirna_id), str(pre['pre_id']))
		    if line not in mat_mirnas_list:
			mat_mirnas_list.append(line)
    """Response"""    
    try:
        return render_template('prediction_analysis/prediction.html',
                               sys_cores = no_system_cores,
                               mat_mirnas_list = mat_mirnas_list)
    except TemplateNotFound:
        return jsonify({ 'error_msg': 'Not Found' }), 404


@prediction_analysis.route('/prediction_analysis_details', methods=['POST'])
def prediction_analysis_details():
    
    if request.method != 'POST':
        return jsonify({ 'error_msg': 'No data available' }), 405
    
    tools = request.form.getlist('tools[]')
    consensus = request.form.getlist('consensus')[0]
    cores = request.form.getlist('cores')[0]
    fasta = request.form.getlist('fa')[0]
    tools_combination = request.form.getlist('combination')[0]
    
    # ==========================================================================
    # Input filtering 
    # ==========================================================================
    tools_f = []
    tools_combination_f = []
    
    for el in tools:
        f_el = re.sub(r'[^0-4]', '', el)
        if int(f_el) in [0,1,2,3,4]:
            tools_f.append(int(f_el))
    
    if tools_combination != '':
	for el in tools_combination.split(','):
	    f_el = re.sub(r'[^0-4]', '', el)
	    if int(f_el) in tools_f:
		tools_combination_f.append(int(f_el))
    
    consensus_f = int(re.sub(r'[^1-5]', '', consensus))
    cores_f = int(re.sub(r'[^0-9]+', '', cores))
    fasta_f = re.sub(r'[^0-9a-z>\.\-\+\|\s\n\r:;]+', '', fasta.replace('<br>','\n'), flags = re.I)
        
    if not tools_f or \
       not isinstance(consensus_f, int) or \
       not isinstance(cores_f, int) or \
       fasta_f == '':
        return jsonify({ 'error_msg': 'No data available' }), 405
    
    mirnas_list = fasta_f.split('\n')
        
    # miRNA object
    mir_obj = iso.MiRNAsProcessing(res.PRE_MAT_MIRNA_PATH,
                                   res.MAT_PRE_ACC_PATH)
    mir_obj.getMiRNAs(mirnas_list)
    
    with open(res.MIRNA_PATH, 'w') as f:
        f.write(json.dumps(mir_obj.mirnas_list, indent=4))
	
    if not os.path.isfile(res.MIRNA_PATH):
        return jsonify({ 'error_msg': 'No data available' }), 405
    
    task = {
        'no_cores': cores_f,
        'consensus': consensus_f,
        'tools': tools_f,
        'tools_combination': tools_combination_f
    }
    
    with open(res.TASK_PATH, 'w') as f:
        f.write(json.dumps(task, indent=4))
    
    if not os.path.isfile(res.TASK_PATH):
        return jsonify({ 'error_msg': 'No data available' }), 405
    
    single_nt_table_content = []
    isomir_table_content = []
    wild_type_table_content = []
    
    errors = []
        
    for mirna in mir_obj.mirnas_list:
        if mirna['variations']['isomir']:
            isomir_table_content.append(mirna)
        if mirna['variations']['substitution']:
            single_nt_table_content.append(mirna)
    
	# No variation found
	if not mirna['variations']['substitution'] and \
	   not mirna['variations']['isomir']:
	    wild_type_table_content.append(mirna)
	    
	header = mirna['header']
	# Get all errors
	if mirna['variations']['error']:
	    for e in mirna['variations']['error']:
		header = header.replace(e['variation'], '<span class="badge badge-pill badge-danger">{0}</span>'.format(e['variation']))
		
	    errors.append({
	        'mirna_id': mirna['identifier'],
	        'errors': mirna['variations']['error'],
	        'header': header
	    })
	    
    """Response"""    
    try:
        return render_template('prediction_analysis/prediction_analysis.html',
                               mirnas_list = mir_obj.mirnas_list,
                               isomir_table_content = isomir_table_content,
                               single_nt_table_content = single_nt_table_content,
	                       wild_type_table_content = wild_type_table_content,
	                       m = mirnas_list,
	                       errors = errors)
    except TemplateNotFound:    
        return jsonify({ 'error_msg': 'Not Found' }), 404


@prediction_analysis.route('/prediction_analysis_sugg_mirna', methods=['POST'])
def prediction_analysis_sugg_mirna():
    
    if request.method != 'POST':
        return jsonify({ 'error_msg': 'No data available' }), 405
        
    data = request.form.to_dict()
    
    if 'mirna' not in data:
        return jsonify({ 'error_msg': 'No data available' }), 405 
    
    # miRNA and pre IDs
    pre_mat_ids = re.sub(r'[^0-9a-zA-Z\.\-_\s:\[\]]+', '', data['mirna']) if data['mirna'] else ''    
    
    match_obj = re.match(r'^(hsa[^\s]+)\s+\[pre\-miR:\s+(hsa[^\s]+)\]$', pre_mat_ids, re.M|re.I)  
    if not match_obj:
	return jsonify({ 'error_msg': 'No data available' }), 405 
    
    mirna_id = match_obj.group(1).strip()
    pre_id = match_obj.group(2).strip()
    
    # Mature miRNAs (Mirbase v20)
    mat_mirnas_list = {}
    if os.path.isfile(res.PRE_MAT_MIRNA_PATH):
        with open(res.PRE_MAT_MIRNA_PATH, 'r') as f:
            mat_mirnas_list = json.load(f) 
            
    if mirna_id not in mat_mirnas_list:
        return jsonify({ 'error_msg': 'No data available' }), 405
    
    mat_seq = ''
    pre_seq = ''
    for pre in mat_mirnas_list[mirna_id]:
	if pre['pre_id'] == pre_id:
	    pre_seq = pre['pre_seq']
	    mat_seq = pre['mature_seq']
	    break
	
    """Response"""    
    try:
        return render_template('prediction_analysis/add_mirna.html',
                               mirna_id = mirna_id,
	                       pre_id = pre_id,
                               mat_seq = mat_seq,
	                       pre_seq = pre_seq)
    except TemplateNotFound:
        return jsonify({ 'error_msg': 'Not Found' }), 404


@prediction_analysis.route('/prediction_analysis_processing_status', methods=['POST'])
def prediction_analysis_processing_status():
    
    if request.method != 'POST':
        return jsonify({ 'error_msg': 'No data available' }), 405
    
    if os.path.isfile(res.JOB_PATH) or \
       os.path.isfile(res.JOB_ENRICH_PATH) or \
       os.path.isfile(res.JOB_ENRICH_COMP_PATH):
	return jsonify({ 'error_msg': 'There is already a running job.' }), 405      
  
    if not os.path.isfile(res.MIRNA_PATH):
        return jsonify({ 'error_msg': 'No data available' }), 405
    
    if not os.path.isfile(res.TASK_PATH):
        return jsonify({ 'error_msg': 'No data available' }), 405    
    
    """Response"""    
    try:
        return render_template('prediction_analysis/prediction_processing_status.html')
    except TemplateNotFound:
        return jsonify({ 'error_msg': 'Not Found' }), 404
    
    
@prediction_analysis.route('/prediction_analysis_status_check', methods=['POST'])
def prediction_analysis_status_check():
    
    if request.method != 'POST':
        return jsonify({ 'error_msg': 'No data available.' }), 405
  
    data = request.form.to_dict()
    
    pid = -1
    if 'pid' in data:
	pid = int(re.sub(r'[^0-9]+', '', data['pid'])) if data['pid'] else -1 
    
    if os.path.isfile(res.JOB_PATH):
	with open(res.JOB_PATH, 'r') as f:
	    job = json.load(f)
	    pid = job['job_pid']
    try:
	if pid == -1:
	    # Check the job
	    if os.path.isfile(res.JOB_PATH):
		with open(res.JOB_PATH, 'r') as f:
		    job = json.load(f)
		    job['job_status'] = 'done'
		    # Update
		    with open(res.JOB_PATH, 'w') as f:
			f.write(json.dumps(job, indent=4))	
	    return jsonify({ 'job_status': 'done' }), 200 
	
	proc = psutil.Process(pid)
	if proc.status() == psutil.STATUS_ZOMBIE:
	    # Check the job
	    if os.path.isfile(res.JOB_PATH):
		with open(res.JOB_PATH, 'r') as f:
		    job = json.load(f)
		    job['job_status'] = 'done'
		    job['job_pid'] = -1
		    # Update
		    with open(res.JOB_PATH, 'w') as f:
			f.write(json.dumps(job, indent=4))	
	    return jsonify({ 'job_status': 'done' }), 200 
	
	return jsonify({ 'job_status': 'running' }), 200 	
    except psutil.NoSuchProcess:
	pass
    
    # Check the job
    if os.path.isfile(res.JOB_PATH):
	with open(res.JOB_PATH, 'r') as f:
	    job = json.load(f)
	    job['job_status'] = 'done'
	    # Update
	    with open(res.JOB_PATH, 'w') as f:
		f.write(json.dumps(job, indent=4))	
    return jsonify({ 'job_status': 'done' }), 200 	
       

@prediction_analysis.route('/exec_prediction_analysis', methods=['POST'])
def exec_prediction_analysis():
    
    if request.method != 'POST':
        return jsonify({ 'error_msg': 'No data available.' }), 405  
    
    if not os.path.isfile(res.MIRNA_PATH):
	return jsonify({ 'error_msg': 'No data available.' }), 405
    
    if not os.path.isfile(res.TASK_PATH):
        return jsonify({ 'error_msg': 'No data available.' }), 405
    
    task = None
    with open(res.TASK_PATH, 'r') as f:
	task = json.load(f)
    
    if not task:
	return jsonify({ 'error_msg': 'No data available.' }), 405
    
    data = request.form.to_dict()
    min_cons = -1    
    if 'min_cons' in data:
	min_cons = int(re.sub(r'[^0-9]+', '', data['min_cons'])) if data['min_cons'] else -1 
    
    if min_cons == -1:
	return jsonify({ 'error_msg': 'No data available.' }), 405
    
    time.sleep(0.1)
    
    command = "/usr/bin/python /app/iso/run.py --cores {0} --tools {1} --comb '[AND={2}]' --pred yes --cons {3}".format(
        task['no_cores'], 
        ','.join(map(str, task['tools'])).replace('"', ''), 
        ','.join(map(str, task['tools_combination'])).replace('"', ''), 
        min_cons)
        
    try:
	task = celery_run_command.delay(command.split()) 
	task_id = task.id
	task_status = celery_run_command.AsyncResult(task_id)
	task_state = task_status.state
	job_details = task.get()
	
	job = {
	    'task_id': task_id, 
	    'job_pid': job_details['pid'],
	    'job_status': job_details['status'],
	    'job_msg': job_details['msg'],
	    'command': command
	}
	
	with open(res.JOB_PATH, 'w') as f:
	    f.write(json.dumps(job, indent=4))	
	
	return jsonify(job), 202

    except ValidationError as e:
	return jsonify({ 'error_msg': e.messages }), 405
    
    
@prediction_analysis.route('/prediction_analysis_cleaner', methods=['POST'])
def prediction_analysis_cleaner():
    
    if request.method != 'POST':
        return jsonify({ 'error_msg': 'No data available' }), 405
		
    remove_file(res.JOB_PATH)
    
    # Get the number of CPU cores
    no_system_cores = psutil.cpu_count()    
    
    pid = request.form.getlist('pid')[0]
    pid = int(re.sub(r'[^0-9]+', '', pid)) if pid else -1  
    
    if pid != -1:
	process = subprocess.Popen(['bash','-c', '/opt/kill_zombies.sh'], stdout=subprocess.PIPE)
	output, error = process.communicate()
	
    pred_path = '{0}/prediction/prediction.json'.format(res.OUTPUT_DIR)
    pred_data = []
    if not os.path.isfile(pred_path):
	return jsonify({ 'error_msg': 'No data available' }), 405

    with open(pred_path, 'r') as f:
	pred_data = json.load(f)	
	    
    show_isomir_table = False
    show_nt_variation_table = False
    show_wt_table = False
        
    for mirna,details in pred_data.iteritems():
	if len(details['variations']['isomir']) > 0:
	    show_isomir_table = True
	if len(details['variations']['substitution']) > 0:
	    show_nt_variation_table = True
	if len(details['variations']['substitution']) == 0 and \
	   len(details['variations']['isomir']) == 0:
	    show_wt_table = True
	    
    zip_path = '{0}/prediction.zip'.format(res.OUTPUT_DIR)
    zip_file = False
    if os.path.isfile(zip_path):
	zip_file = True
	   
    """Response"""    
    try:
        return render_template('prediction_analysis/pred_results.html',
                               pred_data = pred_data,
	                       sys_cores = no_system_cores,
	                       show_isomir_table = show_isomir_table,
	                       show_nt_variation_table = show_nt_variation_table,
	                       show_wt_table = show_wt_table,
	                       zip_file = zip_file)
    except TemplateNotFound:
        return jsonify({ 'error_msg': 'Not Found' }), 404
    

@prediction_analysis.route('/prediction_analysis_kill_all', methods=['POST'])
def prediction_analysis_kill_all():
    
    if request.method != 'POST':
        return jsonify({ 'error_msg': 'No data available.' }), 405   
    
    remove_file(res.JOB_PATH)
    
    pid = request.form.getlist('pid')[0]
    pid = int(re.sub(r'[^0-9]+', '', pid)) if pid else -1  
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
        
    return jsonify({ 'status': 'done' }), 200


@prediction_analysis.route('/prediction_analysis_variation_statistics', methods=['POST'])
def prediction_analysis_variation_statistics():
    
    if request.method != 'POST':
        return jsonify({ 'error_msg': 'No data available' }), 405
    
    data = request.form.to_dict()
    if 'mirna_id' not in data or \
       'vtype' not in data:
	return jsonify({ 'error_msg': 'No data available' }), 405
    
    mirna_id = re.sub(r'[^0-9a-zA-Z\-]+', '', request.form.getlist('mirna_id')[0]) if request.form.getlist('mirna_id')[0] else '' 
    vtype = re.sub(r'[^\sa-zA-Z]+', '', request.form.getlist('vtype')[0]) if request.form.getlist('vtype')[0] else ''
    
    v3p = re.sub(r'[^\-\+0-9]+', '', request.form.getlist('v3p')[0]) if request.form.getlist('v3p')[0] else ''
    v5p = re.sub(r'[^\-\+0-9]+', '', request.form.getlist('v5p')[0]) if request.form.getlist('v5p')[0] else ''
    
    vpos = re.sub(r'[^0-9]+', '', request.form.getlist('vpos')[0]) if request.form.getlist('vpos')[0] else ''
    vref = re.sub(r'[^AGCUT\-\+]+', '', request.form.getlist('vref')[0]) if request.form.getlist('vref')[0] else ''
    vmod = re.sub(r'[^AGCUT\-\+]+', '', request.form.getlist('vmod')[0]) if request.form.getlist('vmod')[0] else ''
    
    if mirna_id == '' or \
       vtype == '' or \
       vtype not in ['Wild Type', 'isomiR', 'Substitution']:
	return jsonify({ 'error_msg': 'No data available' }), 405    
    
    pred_data = {}
    if not os.path.isfile(res.PRE_JSON_PATH):
	return jsonify({ 'error_msg': 'No data available' }), 405
    
    with open(res.PRE_JSON_PATH, 'r') as f:
	pred_data = json.load(f)	
	
    if not pred_data:
	return jsonify({ 'error_msg': 'No data available' }), 405
    
    if mirna_id not in pred_data:
	return jsonify({ 'error_msg': 'No data available' }), 405
    
    stats = pred_data[mirna_id]
    venn_data = []
    table_data = []
    main_data = []
    main_table_data = []
    tool_tab = ''
    
    filename = ''
    
    if 'tool_tab' in data:
	tool_tab = re.sub(r'[^a-zA-Z]+', '', request.form.getlist('tool_tab')[0]) if request.form.getlist('tool_tab')[0] else ''
	tool_tab = tool_tab.replace('predrow', '')
	
	if tool_tab not in ['targets', 'miRanda', 'miRmap', 'TargetScan', 'RNAhybrid', 'PITA']:
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	filename = '{0}_'.format(tool_tab)
    
    if vtype == 'Substitution':
	if int(vpos) == -1 or vref == '' or vmod == '':
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	for s in stats['variations']['substitution']:
	    if s['position'] == int(vpos) and \
	       s['ref'] == vref and \
	       s['mod'] == vmod:
		filename = '{0}_nt_var_{1}_{2}_{3}'.format(mirna_id, 
		                                           s['position'],
		                                           s['ref'],
		                                           s['mod'])
		
		for tool, data in s['predictions'].iteritems():
		    if len(data) > 0:
			if tool != 'targets':
			    venn_data.append({
				'name': tool,
				'data': data
			    })
			else:   
			    main_data = data
			    
		for tool, data in s['support'].iteritems():
		    if len(data) > 0:
			if tool != 'targets':
			    table_data.append({
				'name': tool,
				'data': data
			    })
			else:   
			    main_table_data = data			
    
    elif vtype == 'isomiR':
	if v3p == '' or v5p == '':
	    return jsonify({ 'error_msg': 'No data available' }), 405
	
	for i in stats['variations']['isomir']:
	    if i['5p'] == v5p and i['3p'] == v3p:
		filename = '{0}_isomir_{1}_{2}'.format(mirna_id, 
		                                       i['5p'],
		                                       i['3p'])	
		
		for tool, data in i['predictions'].iteritems():
		    if len(data) > 0:
			if tool != 'targets':
			    venn_data.append({
				'name': tool,
				'data': data
			    })
			else:   
			    main_data = data
			    
		for tool, data in i['support'].iteritems():
		    if len(data) > 0:
			if tool != 'targets':
			    table_data.append({
				'name': tool,
				'data': data
			    })
			else:   
			    main_table_data = data	
    
    elif vtype == 'Wild Type':
	filename = '{0}_wild_type'.format(mirna_id)	
	for tool, data in stats['predictions'].iteritems():
	    if len(data) > 0:
		if tool != 'targets':
		    venn_data.append({
	                'name': tool,
	                'data': data
	            })
		else:   
		    main_data = data
			    
	for tool, data in stats['support'].iteritems():
	    if len(data) > 0:
		if tool != 'targets':
		    table_data.append({
	                'name': tool,
	                'data': data
	            })
		else:   
		    main_table_data = data
		    
    if tool_tab != '':
	filename += '_{0}'.format(tool_tab)
    
    """"""
    """Response"""    
    try:
        return render_template('prediction_analysis/var_statistics.html',
	                       main_data = main_data,
	                       main_table_data =  main_table_data,
                               table_data = table_data,
	                       venn_data = json.dumps(venn_data),
	                       tool_tab = tool_tab,
	                       filename = filename)
    except TemplateNotFound:
        return jsonify({ 'error_msg': 'Not Found' }), 404    
    
    
    
@prediction_analysis.route('/reload_prediction_analysis_details', methods=['POST'])
def reload_prediction_analysis_details():
    
    if request.method != 'POST':
        return jsonify({ 'error_msg': 'No data available' }), 405
    
    if not os.path.isfile(res.MIRNA_PATH):
        return jsonify({ 'error_msg': 'No data available' }), 405
    
    if not os.path.isfile(res.TASK_PATH):
        return jsonify({ 'error_msg': 'No data available' }), 405
    
    # Get the number of CPU cores
    #no_system_cores = psutil.cpu_count()
    
    mirnas_list = {}
    task = {}
    
    with open(res.MIRNA_PATH, 'r') as f:
	mirnas_list = json.load(f)    
       
    single_nt_table_content = []
    isomir_table_content = []
    wild_type_table_content = []
        
    for mirna in mirnas_list:
        if mirna['variations']['isomir']:
            isomir_table_content.append(mirna)
        if mirna['variations']['substitution']:
            single_nt_table_content.append(mirna)
    
	# No variation found
	if not mirna['variations']['substitution'] and \
	   not mirna['variations']['isomir']:
	    wild_type_table_content.append(mirna)
	    
    """Response"""    
    try:
        return render_template('prediction_analysis/prediction_analysis.html',
                               mirnas_list = mirnas_list,
                               isomir_table_content = isomir_table_content,
                               single_nt_table_content = single_nt_table_content,
	                       wild_type_table_content = wild_type_table_content)
    except TemplateNotFound:    
        return jsonify({ 'error_msg': 'Not Found' }), 404
    
    
    
@prediction_analysis.route('/get_variations_db', methods=['POST'])
def get_variations_db():
    
    if request.method != 'POST':
        return jsonify({ 'error_msg': 'No data available' }), 405

    data = request.form.to_dict()
    if 'target' not in data:
	return jsonify({ 'error_msg': 'No data available' }), 405
    
    target = re.sub(r'[^a-z_]+', '', request.form.getlist('target')[0]) if request.form.getlist('target')[0] else ''
    
    if target not in ['edit_db', 'snp_db']:
	return jsonify({ 'error_msg': 'No data available' }), 405
    
    f_path = res.EDITING_DB_JSON_PATH
    if target == 'snp_db':
	f_path = res.SNP_DB_JSON_PATH
	
    if not os.path.isfile(f_path):
	return jsonify({ 'error_msg': 'No data available' }), 405
    
    variations_list = {}
    with open(f_path, 'r') as f:
	variations_list = json.load(f)     
    
    """Response"""    
    try:
        return render_template('prediction_analysis/variations_db.html',
                               variations_list = variations_list,
	                       var_type = target)
    except TemplateNotFound:    
        return jsonify({ 'error_msg': 'Not Found' }), 404
