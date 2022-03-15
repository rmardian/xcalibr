from flask import Blueprint, render_template, request, redirect, url_for, send_from_directory
from werkzeug.utils import secure_filename
import os

model_blueprint = Blueprint('model', __name__)

APP_ROOT = os.path.dirname(os.path.abspath(__file__))
RESOURCES = os.path.join(APP_ROOT, '../resources/')

from app.mod_model.controllers import validateFile, readFile, generateModel, generatePlot, generateReport

@model_blueprint.route('/')
@model_blueprint.route('/index')
@model_blueprint.route('/index.html')
def index():

    return redirect(url_for('index'))

@model_blueprint.route('/dataset')
def dataset():
    return send_from_directory(directory=RESOURCES, path='sample_datasets.zip', as_attachment=True)

@model_blueprint.route('/add_model', methods=['GET', 'POST'])
def add_model():

    if request.method == 'POST':

        if request.form.get('submit')=='Submit':

            with open(os.path.join(RESOURCES, './custom_crn.txt'), 'a') as custom_model:
                custom_model.write(request.form.get('model_name') + '\n')

    return redirect(url_for('index'))

@model_blueprint.route('/run_model', methods=['GET', 'POST'])
def run_model():

    if request.method == 'POST':
        
        fluo = request.files['fluorescence_data']
        od = request.files['od_data']

        if not fluo or not od:
            return "ERROR"
        
        if validateFile(fluo.filename):
            path_to_fluo = os.path.join(RESOURCES, secure_filename(fluo.filename))
            fluo.save(path_to_fluo)

        if validateFile(od.filename):
            path_to_od = os.path.join(RESOURCES, secure_filename(od.filename))
            od.save(path_to_od)

        fluos, ods = readFile(path_to_fluo, path_to_od)
        gates = list(set([i[:-3] for i in fluos.columns.tolist()])) #might need to revisit
        f_dfs, f_sims, f_ts, f_datas = generateModel(fluos, ods, gates)

        with open(os.path.join(RESOURCES, './custom_model.txt'), 'a') as custom_model:
            custom_model.write(request.form.get('model_name') + '\n')

        return redirect(url_for('index'))
    '''
    
    if request.method == 'POST':
        
        fluo = 'marionette_fluo.csv'
        od = 'marionette_od.csv'

        path_to_fluo = os.path.join(RESOURCES, secure_filename(fluo))
        path_to_od = os.path.join(RESOURCES, secure_filename(od))

        fluos, ods = readFile(path_to_fluo, path_to_od)
        gates = list(set([i[:-3] for i in fluos.columns.tolist()])) #might need to revisit
        f_dfs, f_sims, f_ts, f_datas = generateModel(fluos, ods, gates)
        generatePlot(f_ts, f_datas, f_dfs, f_sims, gates)
        generateReport(f_dfs, gates)

        return render_template('model.html', message="Task in queue. You will receive an email for your report when it is completed!")
    '''
    return redirect(url_for('index'))