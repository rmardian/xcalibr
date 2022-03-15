import httplib2
import json

http = httplib2.Http()
url = "http://localhost:9000/api"

def log_in(user, password):

	data = {'username': user,
			'password': password,
			'application': 'xCalibr'
			}
	response, content = http.request(url + "/login", 'POST', json.dumps(data), headers={'Content-Type': 'application/json'})

	session_user = None
	authHeader = None
	if response.status == 200:
		session_user = json.loads(content.decode("utf-8"))['user']['username']
		authHeader = json.loads(content.decode("utf-8"))['authHeader']

	return session_user, authHeader

def create_user(user, password, email, name):

	data = {'username': user,
			'password': password,
			'email': email,
			'name': name,
			'application': 'xCalibr'
			}
	response, content = http.request(url + "/signup", 'POST', json.dumps(data), headers={'Content-Type': 'application/json'})

	session_user = None
	authHeader = None
	if response.status == 200:
		session_user = json.loads(content.decode("utf-8"))['user']['username']
		authHeader = json.loads(content.decode("utf-8"))['authHeader']

	return session_user, authHeader