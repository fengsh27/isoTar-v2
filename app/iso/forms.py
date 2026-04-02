# -*- coding: UTF-8 -*-

from flask_wtf import FlaskForm
from wtforms.validators import InputRequired
from wtforms import StringField, BooleanField

class MirnaAnalysisSuggestionForm(FlaskForm):
    mir = StringField(u'mir', 
                      validators = [InputRequired("Please type a miRNA.")])
