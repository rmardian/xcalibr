from flask import Blueprint, render_template, request, redirect, url_for, send_from_directory, session
from werkzeug.utils import secure_filename
import os

assembly_blueprint = Blueprint('assembly', __name__)

APP_ROOT = os.path.dirname(os.path.abspath(__file__))
RESOURCES = os.path.join(APP_ROOT, '../resources/')

@assembly_blueprint.route('/output', methods=['GET', 'POST'])
def download_output():

    if request.method=='POST':

        opt = request.form.get('assembly')
        if opt=='Parts Setup (Level 0)':

            return send_from_directory(directory=RESOURCES, path='assembly-plan.zip', as_attachment=True)

        return send_from_directory(directory=RESOURCES, path='combinatorial-assembly-4-inputs.zip', as_attachment=True)
