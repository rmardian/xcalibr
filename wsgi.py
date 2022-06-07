from xcalibr import app

if __name__ == '__main__':
    
    '''dev-only'''
    #app.run(host='localhost', port=7777, debug=True, threaded=True)
    
    '''production'''
    #app.run(host='0.0.0.0', port=7777, debug=True)
    app.run(host='0.0.0.0')