from flask import Blueprint, render_template, request, redirect, url_for, send_from_directory
from werkzeug.utils import secure_filename
import os

primer_blueprint = Blueprint('primer', __name__)

APP_ROOT = os.path.dirname(os.path.abspath(__file__))
RESOURCES = os.path.join(APP_ROOT, '../resources/')

@primer_blueprint.route('/output', methods=['GET', 'POST'])
def download_output():

    return send_from_directory(directory=RESOURCES, path='assembly_plan.zip', as_attachment=True)
