from app import app

if __name__ == '__main__':
    
    '''dev-only'''
    app.run(host='localhost', port=5000, debug=True, threaded=True)
    
    '''production'''
    #app.run(host='0.0.0.0')