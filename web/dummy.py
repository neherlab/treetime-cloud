
from flask import Flask, Response, abort, request,render_template,make_response, redirect, url_for, session, send_from_directory, jsonify
from werkzeug import secure_filename
import threading
import numpy as np
import os, random, subprocess, json
app = Flask(__name__)
app.threads = {};
app.debug=True
ALLOWED_EXTENSIONS = ['fasta', 'nwk', 'csv', 'png', 'jpg']
# these modules needed to probe the input data in-place
from Bio import Phylo
import os,sys

dn = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(dn, "static/py"))
import main as main
#
#
def make_id():
    return "".join([chr(random.randint(65,90)) for ii in range(12)])
    #return "HILWQ89P23OI566UHLUIL"

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'GET':
        #return render_template('flask.html', UserId="Dummy")
	return redirect('/' + make_id())

@app.route('/<userid>', methods=['GET', 'POST'])
def index_session(userid):
    if request.method == 'GET':
        return render_template('flask.html', UserId=userid)
    elif request.method == 'POST':
        # save settings
        root = os.path.join("./sessions", userid)
        if not os.path.exists(root):
            os.makedirs(root)

        if 'settings' in request.get_json():
            ss = request.get_json()['settings']
            with open(os.path.join(root, "settings.json"), 'w') as of:
                json.dump(ss, of, True)
        else:
            # error
            return jsonify({'res':'error',
                'message': 'Client-server error: server cannot find proper '
                            'settings in the request'})
            

        # run the function to initialize the steps
        ostate = os.path.join(root, "state.json") # output state file
        steps = main.load_settings(root)
        state = main.write_initial_state(root, steps, ostate) # initial state
        # subprocess: run main script
        app.threads[userid] = threading.Thread(target = main.process, args = (root, steps, state, ostate))
        app.threads[userid].start()
        app.wait_time[userid] = 1

        return jsonify({'res':'OK', 'message':'You can redirect to the wait page'})
        #return redirect('/'+userid+'/progress')
    else:
        pass

@app.route('/<userid>/progress', methods=['GET', 'POST'])
def progress(userid):
    if request.method =='GET':
        return render_template('progress.html', UserId=userid)

@app.route('/<userid>/session_state', methods=['GET', 'POST'])
def get_session_state(userid):
    
    root = os.path.join ('./sessions', userid)
    inf = os.path.join(root, "state.json")
    if not os.path.exists(inf):
        abort(404)
    with open (inf, 'r') as infile:
        json_data = json.load(infile)
        print (json_data)
    #return Response(json.dumps(json_data),  mimetype='application/json')
    return jsonify(**{"steps": json_data})
    
@app.route("/upload/<userid>/file", methods=['GET', 'POST'])
def upload(userid):
    
    folder = os.path.join("./sessions", userid)
    if not os.path.exists(folder):
        os.makedirs(folder)

    if request.method == 'POST':
        
        if 'treefile' in request.files:
            treefile = request.files['treefile']
            treefile.save(os.path.join(folder, "in_tree.nwk"))
        if 'alnfile' in request.files:
            alnfile = request.files['alnfile']
            alnfile.save(os.path.join(folder, "in_aln.fasta"))
        if 'metafile' in request.files:
            metafile = request.files['metafile']
            metafile.save(os.path.join(folder, "in_meta.csv"))
        
        res = {"UploadFile":"OK"}
        
        return jsonify(**res)
            
    #return render_template('results.html', username=username)
    #return "Hello World!"

@app.route('/<userid>/results', methods=['GET', 'POST'])
def results(userid):
    print (userid)
    return render_template('results.html', UserId=userid)

@app.route('/sessions/<userid>/<filename>', methods=['GET', 'POST'])
def send_file(userid, filename):
    uploads = os.path.join(dn, "sessions/" + userid)
    return send_from_directory(uploads, filename) #with open(os.path.join(uploads, filename), 'r') as inf:
    #    json_data = json.load(inf)
    #print (json_data)
    #return jsonify(**json_data)
    

if __name__ == "__main__":
    app.run(port=3000)
