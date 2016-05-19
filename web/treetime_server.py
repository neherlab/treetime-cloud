from flask import Flask, Response, abort, request,render_template,make_response, redirect, url_for, session, send_from_directory, jsonify
from werkzeug import secure_filename
import threading
import numpy as np
import os, random, subprocess, json, shutil
app = Flask(__name__)
app.threads = {};
app.debug=True
ALLOWED_EXTENSIONS = ['fasta', 'nwk', 'csv', 'png', 'jpg']
# these modules needed to probe the input data in-place
from Bio import Phylo
import os,sys
import StringIO

dn = os.path.dirname(os.path.abspath(__file__))
sessions_root = os.path.join(dn , 'sessions')
sys.path.append(os.path.join(dn, "static/py"))

import main as main
from tree_time_config import config_dic as TREETIME_DEFAULT_CONFIG
from tree_time_process import process as TREETIME_PROCESS

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
        cfg = StringIO.StringIO()
        json.dump(TREETIME_DEFAULT_CONFIG, cfg)
        return render_template('flask.html', UserId=userid, Config=TREETIME_DEFAULT_CONFIG)
    elif request.method == 'POST':
        # save settings
        root = os.path.join(sessions_root, userid)
        if not os.path.exists(root):
            os.makedirs(root)
        
        if 'config' in request.get_json():
            ss = request.get_json()['config']
            with open(os.path.join(root, "config.json"), 'w') as of:
                json.dump(ss, of, True)

            app.threads[userid] = threading.Thread(target = TREETIME_PROCESS, args = (root, ss))
            app.threads[userid].start()
            #app.wait_time[userid] = 1
            return jsonify({'res':'OK', 'message':'You can redirect to the wait page'})

        else:
            # error
            return jsonify({'res':'error',
                'message': 'Client-server error: server cannot find proper '
                            'config in the request'})
        
        #return redirect('/'+userid+'/progress')
    else:
        pass

@app.route('/<userid>/example', methods=['GET', 'POST'])
def run_example(userid):
    
    def copy_files(name, root):
        res = {}
        examples = os.path.join(os.path.join(dn, "examples") , name)
        treefile = name + ".nwk"
        alnfile = name + ".fasta"
        metafile = name + ".csv"
        shutil.copyfile(os.path.join(examples, treefile), os.path.join(root, "in_tree.nwk"))
        res["TreeFile"] = treefile
        shutil.copyfile(os.path.join(examples, alnfile), os.path.join(root, "in_aln.fasta"))
        res["AlnFile"] = alnfile
        shutil.copyfile(os.path.join(examples, metafile), os.path.join(root, "in_meta.csv"))
        res["MetaFile"] = metafile
        return res

    if request.method != 'POST':
        abort(404)
    
    root = os.path.join(sessions_root, userid)
    if not os.path.exists(root):
            os.makedirs(root)
    req_data = request.get_json()
    if 'example' not in req_data:
        abort(404)    
    
    res = {}
    
    if req_data['example'] == 'H3N2_NA_20':
        name = 'H3N2_NA_20'
        res = copy_files(name, root)
    elif req_data['example'] == 'H3N2_NA_500':
        name = 'H3N2_NA_500'
        res = copy_files(name, root)
    else:
        abort(404)

    res["UploadFile"] = "OK"
    return  jsonify(**res)

@app.route('/<userid>/progress', methods=['GET', 'POST'])
def progress(userid):
    if request.method =='GET':
        return render_template('progress.html', UserId=userid)

@app.route('/<userid>/session_state', methods=['GET', 'POST'])
def get_session_state(userid):
    
    root = os.path.join (sessions_root, userid)
    inf = os.path.join(root, "session_state.json")
    if not os.path.exists(inf):
        abort(404)
    with open (inf, 'r') as infile:
        json_data = json.load(infile)
        print (json_data)
    #return Response(json.dumps(json_data),  mimetype='application/json')
    return jsonify(**{"steps": json_data})
    
@app.route("/upload/<userid>/file", methods=['GET', 'POST'])
def upload(userid):
    
    folder = os.path.join(sessions_root, userid)
    if not os.path.exists(folder):
        os.makedirs(folder)

    if request.method == 'POST':
        res = {}

        if 'treefile' in request.files:
            treefile = request.files['treefile']
            treefile.save(os.path.join(folder, "in_tree.nwk"))
            res['TreeFile'] = treefile.filename
        if 'alnfile' in request.files:
            alnfile = request.files['alnfile']
            alnfile.save(os.path.join(folder, "in_aln.fasta"))
            res['AlnFile'] = alnfile.filename
        if 'metafile' in request.files:
            metafile = request.files['metafile']
            metafile.save(os.path.join(folder, "in_meta.csv"))
            res['MetaFile'] = metafile.filename
        
        res["UploadFile"] = "OK"
        
        return jsonify(**res)
            
    #return render_template('results.html', username=username)
    #return "Hello World!"

@app.route('/<userid>/results', methods=['GET', 'POST'])
def results(userid):
    print (userid)
    return render_template('results.html', UserId=userid)

@app.route('/sessions/<userid>/<filename>', methods=['GET', 'POST'])
def send_file(userid, filename):
    uploads = os.path.join(sessions_root, userid)
    return send_from_directory(uploads, filename) #with open(os.path.join(uploads, filename), 'r') as inf:
    #    json_data = json.load(inf)
    #print (json_data)
    #return jsonify(**json_data)
    
@app.route('/terms.html/')
def send_terms():
    return render_template('terms.html')

if __name__ == "__main__":
    app.wait_time = {};
    app.threads = {};
    app.debug=True
    app.run(port=3000)
