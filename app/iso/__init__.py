#!/usr/bin/env python
# -*- coding: UTF-8 -*-
from flask import (
    Flask, 
    render_template, 
    request, 
    session, 
    abort, 
    current_app,
    jsonify,
    json,
    Response
)
from webassets.loaders import PythonLoader as PythonAssetsLoader
from flask_wtf.csrf import CSRFError, generate_csrf, validate_csrf
from os.path import abspath, join

from iso import assets
from iso.controllers.main import main
from iso.controllers.homepage import homepage
from iso.controllers.prediction_analysis import prediction_analysis  
from iso.controllers.enrichment_analysis import enrichment_analysis 
from celery import Celery

from iso.extensions import (
    cache,
    assets_env,
    csrf
)


def create_app(object_name):
    """
    An flask application factory, as explained here:
    http://flask.pocoo.org/docs/patterns/appfactories/

    Arguments:
        object_name: the python path of the config object,
                     e.g. iso.settings.ProdConfig
    """

    app = Flask(__name__)

    app.jinja_env.add_extension('jinja2.ext.do')
    
    app.config.from_object(object_name)
    app.config['MAX_CONTENT_LENGTH'] = 4 * 1024 * 1024 * 1024 # 4GB
    
    # initialize cache
    cache.init_app(app) 
    
    celery = Celery('iso.tasks')
    celery.conf.update(app.config)

    # Import and register the different asset bundles
    #assets_env.load_path = abspath(join(app.root_path, '..'))
    #assets_env.register('common_js', common_js)
    assets_env.init_app(app)
    assets_loader = PythonAssetsLoader(assets)
    for name, bundle in assets_loader.load_bundles().items():
        assets_env.register(name, bundle)
        
    csrf.init_app(app)
    
    # register our blueprints
    app.register_blueprint(main)
    app.register_blueprint(homepage)
    app.register_blueprint(prediction_analysis)
    app.register_blueprint(enrichment_analysis)
    
    @app.errorhandler(404)
    def page_not_found(e):
        return Response(json.dumps({ 'error_msg': e.description }), mimetype='application/json'), 404
    
    @app.errorhandler(500)
    def internal_server_error(e):
        return Response(json.dumps({ 'error_msg': e.description }), mimetype='application/json'), 500   
    
    @app.errorhandler(502)
    def bad_gateway(e):
        return Response(json.dumps({ 'error_msg': e.description }), mimetype='application/json'), 502    
    
    @app.errorhandler(CSRFError)
    def csrf_error(e):
        #return render_template('errors/csrf_error.html', reason=e.description), 400
        #return jsonify({ 'error_msg': e.description }), 400 
        return Response(json.dumps({ 'error_msg': e.description }), mimetype='application/json'), 400
   
    return app