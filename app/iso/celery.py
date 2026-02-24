from __future__ import absolute_import
from celery import Celery

'''
celery = Celery('iso',
                broker='amqp://isorabmq:isorabmq@localhost/isorabmq_vhost',
                backend='rpc://',
                include=['iso.tasks'])
'''

celery = Celery('iso',
                broker='amqp://',
                backend='rpc://',
                include=['iso.tasks'])