from flask import Blueprint, render_template, request, redirect, url_for, session, send_from_directory
from werkzeug.utils import secure_filename
import pandas as pd
import numpy as np
import os
import json

design_blueprint = Blueprint('design', __name__)

APP_ROOT = os.path.dirname(os.path.abspath(__file__))
RESOURCES = os.path.join(APP_ROOT, '../resources/')

from app.mod_design.controllers import automateDesign

@design_blueprint.route('/')
@design_blueprint.route('/index')
@design_blueprint.route('/index.html')
def index():

    return redirect(url_for('index'))

@design_blueprint.route('/biodesign', methods=['GET', 'POST'])
def run_design():

    if request.method == 'POST':

        if request.form.get('save')=='Download':

            with open(os.path.join(RESOURCES, './design.json'), 'w') as file:
                file.write(request.form.get('saved_design'))
            return send_from_directory(directory=RESOURCES, path='design.json', as_attachment=True)

        elif request.form.get('upload_circuit')=='Upload':

            with open(os.path.join(RESOURCES, './custom_model.txt'), 'r') as custom_model:
                data = custom_model.readlines()

            with open(os.path.join(RESOURCES, './design.json'), 'r') as saved_design:
                design = saved_design.readlines()

            return render_template('design.html', data=data, saved_design=design)
        
        elif request.form.get('upload_model')=='Upload':
            return redirect(url_for('automate_design'))
        
        elif request.form.get('submit')=='Submit':

            with open(os.path.join(RESOURCES, './model_params.json'), 'r') as json_file:
                model_params = json.load(json_file)

            with open(os.path.join(RESOURCES, './od_params.json'), 'r') as json_file:
                od_params = json.load(json_file)
            
            params_dfs = []
            for p in model_params:
                params_dfs.append(pd.read_json(p))

            od_params_dfs = []
            for p in od_params:
                od_params_dfs.append(pd.read_json(p)['Value'].values)

            return automateDesign(np.arange(0, 1440, 20), params_dfs, od_params_dfs)

    return redirect(url_for('index'))
