# Import flask dependencies
from flask import (
    Blueprint, 
    render_template, 
    jsonify,
    make_response
)
from jinja2 import TemplateNotFound
import requests
import json

homepage = Blueprint('homepage', 
                     __name__, 
                     template_folder = 'templates')


@homepage.route('/')
def index():
    try:
        return render_template('homepage/index.html')
    except TemplateNotFound:
        return render_template('errors/404.html'), 404
