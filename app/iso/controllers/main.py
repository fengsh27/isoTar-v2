from flask import (
    Blueprint, 
    render_template, 
    flash, 
    request, 
    redirect, 
    url_for,
    jsonify
)

from jinja2 import TemplateNotFound
from iso.extensions import cache
from werkzeug import secure_filename
import os
import os.path
from iso.settings import AppResources as res
import subprocess
import shutil
import json
from iso.tasks import celery_run_command
from marshmallow import ValidationError
import psutil


main = Blueprint('main', 
                 __name__, 
                 template_folder = 'templates')

def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in res.ALLOWED_EXTENSIONS


def execute_command_line(command):
    retcode = 0
    try:
	retcode = subprocess.check_output(command, shell = True)
    except OSError as e:
	pass
	
    return retcode 


def remove_directory(path):
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
	
	
def remove_directory_content(path):
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
                        pass        
                else:
                    pass
	    

@main.route('/doc')
@cache.cached(timeout=1000)
def doc():
    try:
        return render_template('main/doc.html')
    except TemplateNotFound:
        return jsonify({ 'error_msg': 'Not Found' }), 404
    
    
@main.route('/about')
def about():
    try:
        return render_template('main/about.html')
    except TemplateNotFound:
        return jsonify({ 'error_msg': 'Not Found' }), 404
    
    
@main.route('/contact')
def contact():
    try:
        return render_template('main/contact.html')
    except TemplateNotFound:
        return jsonify({ 'error_msg': 'Not Found' }), 404
    

@main.route('/editing_db')
@cache.cached(timeout=1000)
def editing_db():
    
    variations_db = []
    if os.path.isfile(res.EDITING_DB_JSON_PATH):
	with open(res.EDITING_DB_JSON_PATH, 'r') as f:
	    variations_db = json.load(f)    
    try:
        return render_template('main/db.html',
	                       variations_db = variations_db,
	                       var_type = 'edit_db',
	                       filename = 'Known A-to-I Mature miRNA Editing Sites')
    except TemplateNotFound:
        return jsonify({ 'error_msg': 'Not Found' }), 404
    

@main.route('/snp_db')
@cache.cached(timeout=1000)
def snp_db():
    
    variations_db = []
    if os.path.isfile(res.SNP_DB_JSON_PATH):
	with open(res.SNP_DB_JSON_PATH, 'r') as f:
	    variations_db = json.load(f)    
    try:
        return render_template('main/db.html',
	                       variations_db = variations_db,
	                       var_type = 'snp_db',
	                       filename = 'SNPs in Mature miRNA')
    except TemplateNotFound:
        return jsonify({ 'error_msg': 'Not Found' }), 404
    
    
    
    
    
@main.route('/work_uploader', methods=['POST'])
def upload_file():
    if request.method != 'POST':
        return jsonify({ 'error_msg': 'No file uploaded' }), 405 
    
    if os.path.isfile(res.JOB_PATH) or \
       os.path.isfile(res.JOB_ENRICH_PATH) or \
       os.path.isfile(res.JOB_ENRICH_COMP_PATH):
	return jsonify({ 'error_msg': 'There is already a running job.' }), 405
    
    files = []    
    # Check if the post request has the file part
    if 'files[]' not in request.files:  
        return jsonify({ 'error_msg': 'No file uploaded' }), 405
    
    file = request.files['files[]']
    
    # Submit a empty part without filename
    if file.filename == '':        
        return jsonify({ 'error_msg': 'Empty file' }), 405
    
    if not file or not allowed_file(file.filename):
        return jsonify({ 'error_msg': 'Bad file extension' }), 405
    
    filename = secure_filename(file.filename)
    file.save(os.path.join(res.UPLOAD_PATH, filename))
    f = {
        'name': filename,
    }
    
    if not os.path.isfile(os.path.join(res.UPLOAD_PATH, filename)):
	return jsonify({ 'error_msg': 'Uploaded file is not valid' }), 405     
    
    # Remove all previous work
    remove_directory_content(res.OUTPUT_DIR)
    
    files.append(f)
    return jsonify(files=files)  


@main.route('/prepare_uploaded_file', methods=['POST'])
def prepare_uploaded_file():
    if request.method != 'POST':
        return jsonify({ 'error_msg': 'No data available' }), 405 
    
    if not os.path.isfile(os.path.join(res.UPLOAD_PATH, 'output.zip')):  
        return jsonify({ 'error_msg': 'No uploaded file found' }), 405
    
    # Unzip the file
    command = '/usr/bin/unzip -o {0} -d {1}'.format(
        os.path.join(res.UPLOAD_PATH, 'output.zip'),
        os.path.join(res.OUTPUT_DIR.replace('output', '')))    
       
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


@main.route('/upload_work_check', methods=['POST'])
def upload_work_check():
    
    pid = -1
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
	    
	    # Remove the zip file
	    if os.path.isfile(os.path.join(res.UPLOAD_PATH, 'output.zip')): 
		remove_directory(os.path.join(res.UPLOAD_PATH, 'output.zip'))
		
	    # Remove the job file
	    if os.path.isfile(os.path.join(res.OUTPUT_DIR, res.JOB_PATH)): 
		remove_directory(os.path.join(res.OUTPUT_DIR, res.JOB_PATH))
	
	    # Remove zombie processes
	    process = subprocess.Popen(['bash','-c', '/opt/kill_zombies.sh'], stdout=subprocess.PIPE)
	    output, error = process.communicate()
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
			
	    # Remove the zip file
	    if os.path.isfile(os.path.join(res.UPLOAD_PATH, 'output.zip')): 
		remove_directory(os.path.join(res.UPLOAD_PATH, 'output.zip'))
		
	    # Remove the job file
	    if os.path.isfile(os.path.join(res.OUTPUT_DIR, res.JOB_PATH)): 
		remove_directory(os.path.join(res.OUTPUT_DIR, res.JOB_PATH))
	
	    # Remove zombie processes
	    process = subprocess.Popen(['bash','-c', '/opt/kill_zombies.sh'], stdout=subprocess.PIPE)
	    output, error = process.communicate()
	    
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
    
    # Remove the zip file
    if os.path.isfile(os.path.join(res.UPLOAD_PATH, 'output.zip')): 
	remove_directory(os.path.join(res.UPLOAD_PATH, 'output.zip'))
    
    # Remove the job file
    if os.path.isfile(os.path.join(res.OUTPUT_DIR, res.JOB_PATH)): 
	remove_directory(os.path.join(res.OUTPUT_DIR, res.JOB_PATH))
	
    # Remove zombie processes
    process = subprocess.Popen(['bash','-c', '/opt/kill_zombies.sh'], stdout=subprocess.PIPE)
    output, error = process.communicate()

    return jsonify({ 'job_status': 'done' }), 200 



@main.route('/work_dowloader')
def download_file():
    
    if not os.path.isdir(res.OUTPUT_DIR):
	return jsonify({ 'error_msg': 'No data available' }), 405  
    
    if os.path.isfile(res.JOB_PATH) or \
       os.path.isfile(res.JOB_ENRICH_PATH) or \
       os.path.isfile(res.JOB_ENRICH_COMP_PATH):
	return jsonify({ 'error_msg': 'There is already a running job.' }), 405
    
    # Remove the previous zip file
    if os.path.isfile(os.path.join(res.OUTPUT_DIR.replace('output', ''), 
                                   'output.zip')):
	remove_directory(os.path.join(res.OUTPUT_DIR.replace('output', ''), 
	                              'output.zip')) 
    
    command = '/usr/bin/zip -r output.zip output'

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
    
    return jsonify({ 'status': 'done' }), 200


@main.route('/download_work_check', methods=['POST'])
def download_work_check():
    
    pid = -1
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
	
	    # Remove the job file
	    if os.path.isfile(os.path.join(res.OUTPUT_DIR, res.JOB_PATH)): 
		remove_directory(os.path.join(res.OUTPUT_DIR, res.JOB_PATH))
	
	    # Remove zombie processes
	    process = subprocess.Popen(['bash','-c', '/opt/kill_zombies.sh'], stdout=subprocess.PIPE)
	    output, error = process.communicate()
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

	    # Remove the job file
	    if os.path.isfile(os.path.join(res.OUTPUT_DIR, res.JOB_PATH)): 
		remove_directory(os.path.join(res.OUTPUT_DIR, res.JOB_PATH))
	
	    # Remove zombie processes
	    process = subprocess.Popen(['bash','-c', '/opt/kill_zombies.sh'], stdout=subprocess.PIPE)
	    output, error = process.communicate()
	    
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
   
    # Remove the job file
    if os.path.isfile(os.path.join(res.OUTPUT_DIR, res.JOB_PATH)): 
	remove_directory(os.path.join(res.OUTPUT_DIR, res.JOB_PATH))
	
    # Remove zombie processes
    process = subprocess.Popen(['bash','-c', '/opt/kill_zombies.sh'], stdout=subprocess.PIPE)
    output, error = process.communicate()

    return jsonify({ 'job_status': 'done' }), 200 
     
