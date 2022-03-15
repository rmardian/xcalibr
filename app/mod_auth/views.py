from flask import Blueprint, render_template, request, redirect, url_for, session

auth_blueprint = Blueprint('auth', __name__)

from app.mod_auth.controllers import log_in, create_user

@auth_blueprint.route('/auth')
def auth():
    return 'Authentication module!'

@auth_blueprint.route("/register", methods=['POST'])
def register():

	if request.method == 'POST':
		session_user, authHeader = create_user(request.form['username'], request.form['password'], request.form['email'], request.form['name'])
		if session_user != None:
			session['user'] = session_user
			session['logged_in'] = True
			session['authHeader'] = authHeader
			return redirect(url_for('search'))
		return render_template('login.html', message="Error: Username or email has been taken!")
	return render_template('login.html')

@auth_blueprint.route("/login.html", methods=['GET', 'POST'])
@auth_blueprint.route("/login", methods=['GET', 'POST'])
def login():

	if request.method == 'POST':
		#send request for authentication
		#session_user, authHeader = log_in(request.form['username'], request.form['password'])
		#if session_user != None:
		#	session['user'] = session_user
		#	session['logged_in'] = True
		#	session['authHeader'] = authHeader
		#	return redirect(url_for('index'))
		#return render_template('login.html', message="Login error! Please check your username and password")
		
		session['user'] = request.form['username']
		session['logged_in'] = True
		return redirect(url_for('index'))
	return render_template('login.html')

@auth_blueprint.route('/logout.html', methods=['GET', 'POST'])
@auth_blueprint.route('/logout', methods=['GET', 'POST'])
def logout():
	session.pop('user', None)
	session.pop('logged_in', None)
	session.pop('authHeader', None)
	return redirect(url_for('auth.login'))