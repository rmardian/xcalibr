from flask import Flask, render_template, request, redirect, url_for, session, send_from_directory
import os

app = Flask(__name__)
app.secret_key = "06251987"

APP_ROOT = os.path.dirname(os.path.abspath(__file__))
RESOURCES = os.path.join(APP_ROOT, './resources/')

from app.mod_auth.views import auth_blueprint
from app.mod_model.views import model_blueprint
from app.mod_design.views import design_blueprint
from app.mod_primer.views import primer_blueprint

app.register_blueprint(auth_blueprint, url_prefix='/auth')
app.register_blueprint(model_blueprint, url_prefix='/generate_model')
app.register_blueprint(design_blueprint, url_prefix='/automate_design')
app.register_blueprint(primer_blueprint, url_prefix='/primer')

@app.route("/")
@app.route("/index.html")
@app.route("/index")
def index():

	#session['user'] = 'xcalibr'
	#session['logged_in'] = True

	if 'user' not in session or not session['logged_in']:
		return redirect(url_for('auth.login'))
	return redirect(url_for('generate_model'))

@app.route('/model.html', methods=['GET', 'POST'])
@app.route('/model', methods=['GET', 'POST'])
def generate_model():

	if 'user' not in session or not session['logged_in']:
		return redirect(url_for('auth.login'))

	with open(os.path.join(RESOURCES, './custom_crn.txt'), 'r') as custom_model:
		data = custom_model.readlines()

	return render_template('model.html', data=data)

@app.route('/design.html', methods=['GET', 'POST'])
@app.route('/design', methods=['GET', 'POST'])
def automate_design():

	if 'user' not in session or not session['logged_in']:
		return redirect(url_for('auth.login'))

	with open(os.path.join(RESOURCES, './custom_model.txt'), 'r') as custom_model:
		data = custom_model.readlines()

	return render_template('design.html', data=data)

@app.route('/simulation.html', methods=['GET', 'POST'])
@app.route('/simulation', methods=['GET', 'POST'])
def simulate_design():

	if 'user' not in session or not session['logged_in']:
		return redirect(url_for('auth.login'))

	with open(os.path.join(RESOURCES, './custom_model.txt'), 'r') as custom_model:
		data = custom_model.readlines()

	return render_template('simulation.html', data=data)

@app.route('/primer.html', methods=['GET', 'POST'])
@app.route('/primer', methods=['GET', 'POST'])
def compose_primer():

	if 'user' not in session or not session['logged_in']:
		return redirect(url_for('auth.login'))

	return render_template('primer.html')

@app.route('/custom.html', methods=['GET', 'POST'])
@app.route('/custom', methods=['GET', 'POST'])
def custom_model():

	if 'user' not in session or not session['logged_in']:
		return redirect(url_for('auth.login'))

	return render_template('custom.html')

@app.errorhandler(404) 
def not_found(e): 
	return render_template("page_404.html") 
