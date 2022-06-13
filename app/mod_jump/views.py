from flask import Blueprint, render_template, request, redirect, url_for, session, send_from_directory, Markup
from werkzeug.utils import secure_filename
import pandas as pd
import numpy as np
import os
import json

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

jump_blueprint = Blueprint('jump', __name__)

APP_ROOT = os.path.dirname(os.path.abspath(__file__))
RESOURCES = os.path.join(APP_ROOT, '../resources/')
INPUT = os.path.join(RESOURCES, './input/')
OUTPUT = os.path.join(RESOURCES, './output/')

from app.mod_jump.controllers import handle_file, domesticate_parts, simulate_assembly, visualize_construct

@jump_blueprint.route('/')
@jump_blueprint.route('/index')
@jump_blueprint.route('/index.html')
def index():

    return redirect(url_for('jump_home'))

@jump_blueprint.route('/domestication', methods=['GET', 'POST'])
def run_part_domestication():

    if request.method == 'POST':

        title = 'Part domestication'
        parts_req = request.files['parts']
        map_req = request.files['map']
        plasmids_req = request.files['plasmids']

        #if either mandatory files not provided
        if not parts_req or not plasmids_req:
            return render_template('jump-report.html', title=title, header="File missing, task not executed.", status='Failed')

        parts_file = handle_file(parts_req, 'tabular')
        map_file = handle_file(map_req, 'tabular') if map_req else None
        plasmids_file = handle_file(plasmids_req, 'sequence')
        
        #if any file extension invalid
        if not parts_file or not plasmids_file or (map_req and not map_file):
            return render_template('jump-report.html', title=title, header="Upload error, task not executed.", status='Failed')

        #if everything checks out
        parts = pd.read_csv(INPUT + '{}'.format(parts_file[0])) #always at index 0
        parts['sequence'] = parts['sequence'].str.rstrip().str.upper()

        mapping = pd.read_csv(INPUT + '{}'.format(map_file[0])) if map_file else None

        fastas = [f for f in plasmids_file if not f.startswith('__')] #just to make sure valid fasta files
        plasmids = pd.DataFrame([(p.id, str(p.seq)) for fasta in fastas \
                                for p in list(SeqIO.parse(INPUT + '{}'.format(fasta), 'fasta'))], \
                                columns=['name', 'sequence'])
        plasmids['sequence'] = plasmids['sequence'].str.rstrip().str.upper()

        message = domesticate_parts(parts, mapping, plasmids)
        return render_template('jump-report.html', title=title, header="Task successfully executed.", message=Markup(message), status='OK')

    return redirect(url_for('jump_home'))

@jump_blueprint.route('/simulation', methods=['GET', 'POST'])
def run_assembly_simulation():

    if request.method == 'POST':

        title = 'Assembly simulation'
        enzyme = request.form.get('enzyme')
        plan_req = request.files['plan']
        map_req = request.files['map']
        plasmids_req = request.files['plasmids']

        #if either mandatory files not provided
        if not plan_req or not plasmids_req:
            return render_template('jump-report.html', title=title, header="File missing, task not executed.", status='Failed')

        plan_file = handle_file(plan_req, 'tabular')
        map_file = handle_file(map_req, 'tabular') if map_req else None
        plasmids_file = handle_file(plasmids_req, 'sequence')

        #if any file extension invalid
        if not plan_file or not plasmids_file or (map_req and not map_file):
            return render_template('jump-report.html', title=title, header="Upload error, task not executed.", status='Failed')

        #if everything checks out
        plan = pd.read_csv(INPUT + '{}'.format(plan_file[0]))

        mapping = pd.read_csv(INPUT + '{}'.format(map_file[0])) if map_file else None
        
        fastas = [f for f in plasmids_file if not f.startswith('__')] #just to make sure valid fasta files
        plasmids = pd.DataFrame([(p.id, str(p.seq), fasta) for fasta in fastas \
                            for p in list(SeqIO.parse(INPUT + '{}'.format(fasta), 'fasta'))], \
                            columns=['name', 'sequence', 'level'])
        plasmids['level'] = plasmids['level'].str.split('.', expand=True)[0]
        plasmids['sequence'] = plasmids['sequence'].str.rstrip().str.upper()

        message = simulate_assembly(plan, mapping, plasmids, enzyme)
        
        #details = "Part {} cannot be amplified, please find correct templates to amplify them from.".format(missing)
        return render_template('jump-report.html', title=title, header="Task successfully executed.", message=Markup(message), status='OK')

    return redirect(url_for('jump_home'))

@jump_blueprint.route('/visualization', methods=['GET', 'POST'])
def run_construct_visualization():

    if request.method == 'POST':

        plan = request.files['plan']
        parts = request.files['parts']

        #if either mandatory files not provided
        if not plan or not parts:
            return render_template('jump-report.html', title='Construct visualization', header="File missing, task not executed.", status='Failed')

        plan_file = handle_file(plan, 'tabular')
        parts_file = handle_file(parts, 'tabular')

        #if any file extension invalid
        if not plan_file or not parts_file:
            return render_template('jump-report.html', title='Construct visualization', header="Upload error, task not executed.", status='Failed')

        message = visualize_construct(plan_file, parts_file)
        
        return render_template('jump-report.html', title='Construct visualization', header="Task successfully executed.", status='OK')

    return redirect(url_for('jump_home'))

@jump_blueprint.route('/output', methods=['GET', 'POST'])
def download_output():

    if request.method=='POST':

        return send_from_directory(directory=OUTPUT, path='output.zip', as_attachment=True)

    return redirect(url_for('jump_home'))
