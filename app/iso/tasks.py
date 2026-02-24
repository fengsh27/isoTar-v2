from __future__ import absolute_import
from iso.celery import celery
import time
import subprocess

@celery.task(ignore_result=False)
def celery_run_command(command):
    #time.sleep(1)
    d = None
    try: 
	cmd = subprocess.Popen(command,
	                       cwd = '/app/iso/static/public')	
	
	d = { 'status' : 'running', 'msg': '', 'pid': cmd.pid }
    except subprocess.CalledProcessError as e: 
	d = { 'status' : 'error', 
	      'msg': 'An error occurred during the analysis processing.',
	      'pid': -1
	}
	pass
    
    return d
    
    