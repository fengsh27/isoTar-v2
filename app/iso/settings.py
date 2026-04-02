import os
import hashlib

class Config(object):
    SECRET_KEY = '76wtfw3vbejfi8gyhrenwksinunrwiuegnuiwxnhiewquy3oih3c32r8thv0v'
    WTF_CSRF_SECRET_KEY = '76wtfw3vbejfi8gyhrenwksinunrwiuegnuiwxnhiewquy3oih3c32r8thv0v'
    WTF_CSRF_TIME_LIMIT = 84600
    TIME_LIMIT = 3600
    #CELERY_BROKER_URL = 'amqp://isorabmq:isorabmq@localhost/isorabmq_vhost'
    CELERY_BROKER_URL = 'amqp://'
    CELERY_RESULT_BACKEND = 'rpc://'
    # thats where celery will store scheduled tasks in case you restart the broker:
    CELERYD_STATE_DB = "/app/iso/celery_status" 
    CELERY_ACCEPT_CONTENT = ['json']
    CELERY_TASK_SERIALIZER = 'json'
    CELERY_RESULT_SERIALIZER = 'json'  
    CELERY_SEND_TASK_SENT_EVENT = True


class ProdConfig(Config):
    ENV = 'prod'
    CACHE_TYPE = 'simple'
    WTF_CSRF_ENABLED = False
    TIME_LIMIT = 3600
    ASSETS_DEBUG = True
    #SESSION_COOKIE_SECURE = False


class DevConfig(Config):
    ENV = 'dev'
    DEBUG = True
    DEBUG_TB_INTERCEPT_REDIRECTS = False
    CACHE_TYPE = 'null'
    ASSETS_DEBUG = True
    WTF_CSRF_ENABLED = False
    TIME_LIMIT = 3600

    
class AppResources():
    MAT_MIRNA_PATH = '/app/iso/static/public/resources/mature_mirna.json'
    PRE_MAT_MIRNA_PATH = '/app/iso/static/public/resources/marute_pre_mirna.json'
    MIRNA_PATH = '/app/iso/static/public/output/mirnas.json'
    TARGET_PAT = '/app/iso/static/public/resources/3utr.fa'
    OUTPUT_DIR = '/app/iso/static/public/output'
    GO_PATH = '/app/iso/static/public/resources/go.json' 
    MAP_PATH = '/app/iso/static/public/resources/targets_go_map.tab'
    TASK_PATH = '/app/iso/static/public/output/task.json'
    JOB_PATH = '/app/iso/static/public/output/job.json'
    JOB_ENRICH_PATH = '/app/iso/static/public/output/job_go.json'
    JOB_ENRICH_COMP_PATH = '/app/iso/static/public/output/job_go_comp.json'
    MAT_PRE_ACC_PATH = '/app/iso/static/public/resources/mat_pre_acc.json'
    NM_TO_GENEID_PATH = '/app/iso/static/public/resources/nm_to_geneid.tab'
    PRE_JSON_PATH = '/app/iso/static/public/output/prediction/prediction.json'    
    ENRICH_JSON_PATH = '/app/iso/static/public/output/enrichment/enrichment.json'
    ENRICH_COMP_JSON_PATH = '/app/iso/static/public/output/enrichment_comp/enrichment.json'
    OBO_DAG_PATH = '/app/iso/static/public/resources/go-basic.obo'
    GENE_2_GO_PATH = '/app/iso/static/public/resources/gene2go'  
    UPLOAD_PATH = '/tmp'
    ALLOWED_EXTENSIONS = set(['zip'])
    EDITING_DB_JSON_PATH = '/app/iso/static/public/resources/editing_db.json'
    SNP_DB_JSON_PATH = '/app/iso/static/public/resources/snp_db.json'
    PATHWAY_PATH = '/app/iso/static/public/resources/pathways.json'
