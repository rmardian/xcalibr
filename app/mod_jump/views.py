from flask import Blueprint, render_template, request, redirect, url_for, session, send_from_directory
from werkzeug.utils import secure_filename
import pandas as pd
import numpy as np
import os
import json

jump_blueprint = Blueprint('jump', __name__)

APP_ROOT = os.path.dirname(os.path.abspath(__file__))
RESOURCES = os.path.join(APP_ROOT, '../resources/')

from app.mod_jump.controllers import handle_file, domesticate_parts, simulate_assembly

@jump_blueprint.route('/')
@jump_blueprint.route('/index')
@jump_blueprint.route('/index.html')
def index():

    return redirect(url_for('jump_home'))

@jump_blueprint.route('/domestication', methods=['GET', 'POST'])
def run_part_domestication():

    if request.method == 'POST':

        parts = request.files['parts']
        maps = request.files['maps']
        plasmids = request.files['plasmids']

        #if either mandatory files not provided
        if not parts or not plasmids:
            return render_template('jump-report.html', title='Part domestication', message="File missing, task not executed.")

        parts_file = handle_file(parts, 'tabular')
        maps_file = None
        if maps:
            maps_file = handle_file(maps, 'tabular')
        plasmids_file = handle_file(plasmids, 'sequence')

        #if any file extension invalid
        if not parts_file or not plasmids_file or (maps and not maps_file):
            return render_template('jump-report.html', title='Part domestication', message="Upload error, task not executed.")

        message = domesticate_parts(parts_file, maps_file, plasmids_file)
        
        #details = "Part {} cannot be amplified, please find correct templates to amplify them from.".format(missing)
        return render_template('jump-report.html', title='Part domestication', message="Task successfully executed.")

    return redirect(url_for('jump_home'))

@jump_blueprint.route('/simulation', methods=['GET', 'POST'])
def run_assembly_simulation():

    if request.method == 'POST':

        enzyme = request.form.get('enzyme')
        plan = request.files['plan']
        maps = request.files['maps']
        plasmids = request.files['plasmids']

        #if either mandatory files not provided
        if not plan or not plasmids:
            return render_template('jump-report.html', title='Assembly simulation', message="File missing, task not executed.")

        plan_file = handle_file(plan, 'tabular')
        maps_file = None
        if maps:
            maps_file = handle_file(maps, 'tabular')
        plasmids_file = handle_file(plasmids, 'sequence')

        #if any file extension invalid
        if not plan_file or not plasmids_file or (maps and not maps_file):
            return render_template('jump-report.html', title='Assembly simulation', message="Upload error, task not executed.")

        message = simulate_assembly(plan_file, maps_file, plasmids_file, enzyme)
        
        #details = "Part {} cannot be amplified, please find correct templates to amplify them from.".format(missing)
        return render_template('jump-report.html', title='Assembly simulation', message="Task successfully executed.")

    return redirect(url_for('jump_home'))
