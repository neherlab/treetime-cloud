from flask import Flask, Response, abort, request,render_template,make_response, redirect, url_for, session, send_from_directory, jsonify
from werkzeug import secure_filename
import threading
import numpy as np
import os, random, subprocess, json, shutil
from Bio import Phylo
import os,sys
import StringIO

app = Flask(__name__)
app.threads = {};
app.debug=True
ALLOWED_EXTENSIONS = ['fasta', 'nwk', 'csv', 'png', 'jpg']

# html theme is controlled by the server
html_theme_css = "http://bootswatch.com/flatly/bootstrap.css"

# set the server directories
dn = os.path.dirname(os.path.abspath(__file__))
sessions_root = os.path.join(dn , 'sessions')

# import additional modules
# NOTE the path is relative to this script, should be appended to sys.path beforehand
sys.path.append(os.path.join(dn, "static/py"))
from tree_time_config import treetime_webconfig, treeanc_webconfig
from tree_time_process import run_treetime as RUN_TREETIME
from tree_time_process import run_treeanc as RUN_TREEANC

def make_id():
    """
    Create new userID, which will be used to identigy session.
    """
    return "".join([chr(random.randint(65,90)) for ii in range(12)])

@app.route('/', methods=['GET', 'POST'])
def index():
    """
    Client request for the welcome page
    """
    if request.method == 'GET':
        return render_template('welcome.html', html_theme=html_theme_css)
    return None

@app.route('/ancestral_reconstruction_request', methods=['GET', 'POST'])
def ancestral_reconstruction_request():
    """
    Response on user press 'Ancestral reconstruction' button. It creates new user
    ID and sends it back hence allowing the client to redirect to the ancestral
    reconstruction welcome page
    """
    userid = "anc_" + make_id();
    print ("Ancestral reconstruction request " + request.method + "  user_id: " + userid)

    data = {
            "UserId": userid,
            "config": {}
        }

    response = app.response_class(
            response=json.dumps(data),
            status=200,
            mimetype='application/json'
    )
    return response


@app.route('/ancestral/<userid>', methods=['GET', 'POST'])
def ancestral_reconstruction_welcome(userid):
    """
    Render ancestral reconstruction welcome page and send it to the client
    """
    if request.method == 'GET':
        return render_template('welcome_ancestral_reconstruction.html', UserId=userid, Config=treeanc_webconfig, html_theme=html_theme_css)
    else:
        pass

@app.route('/treetime_request', methods=['GET', 'POST'])
def treetime_request():
    """
    Server response on the user 'TreeTime reconstruction' button pressed. It creates
    new UserID with 'tt_' prefix and sends it back to the web client hence allowing
    the client to redirect to the treetime configuration page.
    """
    userid = "tt_" + make_id();
    print ("treetime request: " + request.method + "  user_id: " + userid)
    data = {
        "UserId": userid,
        "config": {}
    }
    response = app.response_class(
            response=json.dumps(data),
            status=200,
            mimetype='application/json'
    )
    return response

@app.route('/treetime/<userid>', methods=['GET', 'POST'])
def render_treetime_welcome(userid):
    """
    Process user request to render the treetime welcome page. Render the HTML template
    and send it back to the user.
    """
    if request.method == 'GET':
        return render_template('welcome_treetime.html', UserId=userid, Config=treetime_webconfig, html_theme=html_theme_css)
    else:
        pass

@app.route('/upload/treetime/<userid>/file', methods=['POST'])
def upload_treetime_file(userid):
    """
    Upload user-provided file for treetime reconstruction. Get the file type from
    the POST request metadata. Assigns the unified filename based on the file type
    (alignment, tree or metadata file), and stores it in the 'sessions/<userid>'
    folder.
    """
    # define the destination folder
    folder = os.path.join(sessions_root, userid)
    if not os.path.exists(folder):
        os.makedirs(folder)


    if request.method == 'POST':
        res = {}
        # todo validate the content of uploaded files
        # switch filetype and save files under unified names
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


@app.route('/treetime/<userid>/example', methods=['POST'])
def upload_example(userid):
    """
    User request to run of one of the predefined datasets (examples). Based on
    the example name, extract files from the prepared folder and copy them to
    the 'sessions/<userid>' folder.
    """
    def copy_files(name, root):
        """
        Copy example dataset files to the destination folder.
        Args:
         - name(str): name of the example folder (should repeate the name of the
         example files inside).
         - root(str): destination (sessions/<userid>) folder.

        Returns:
         - res(dic): dictionary of the original filenames for the particular example.

        NOTE:
         - The names of the files must repeat the name of the parent folder and
         should be the same as the name of the example, but have appropriate extensions,
         which define the type of the file. I.e., the server expects the
         examples to be organized as follows:
           /examples
            |
            |
            |-/example1
            |  |
            |  |
            |  |-example1.csv
            |  |-example1.nwk
            |  |-example1.fasta
            |
            |
            |-/example2
            |  |
            |  |
            |  |-example2.csv
            |  |-example2.nwk
            |  |-example2.fasta
            |
            |
            ...

        """
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
    # switch the example name, relate to the example directory, copy files
    if req_data['example'] == 'H3N2_NA_20':
        name = 'H3N2_NA_20'
        res = copy_files(name, root)
    elif req_data['example'] == 'H3N2_NA_500':
        name = 'H3N2_NA_500'
        res = copy_files(name, root)
    elif req_data['example'] == 'zika_65':
        name = 'zika_65'
        res = copy_files(name, root)
    elif req_data['example'] == 'H3N2_HA_100':
        name = 'H3N2_HA_100'
        res = copy_files(name, root)
    elif req_data['example'] == 'HIV_RT_200':
        name = 'HIV_RT_200'
        res = copy_files(name, root)
    elif req_data['example'] == 'HIV_p17_200':
        name = 'HIV_p17_200'
        res = copy_files(name, root)
    elif req_data['example'] == 'ebola':
        name = 'ebola'
        res = copy_files(name, root)
    else:
        abort(404)

    res["UploadFile"] = "OK"
    return  jsonify(**res)

@app.route('/treetime/<userid>/run', methods=['POST'])
def run_treetime(userid):
    """
    The server response on user press 'RuntreeTime' buton. Gets the actual treetime
    configuration from the post request metadata, saves it to the 'sessions/<userid>'
    forlder, starts a new thread to run treetime computations.
    """
    if (request.method != "POST"):
        abort(404)

    # save settings
    root = os.path.join(sessions_root, userid)
    if not os.path.exists(root):
        os.makedirs(root)

    #save config, run treetime in a separate thread
    if 'config' in request.get_json():
        ss = request.get_json()['config']

        with open(os.path.join(root, "config.json"), 'w') as of:
            json.dump(ss, of, True)
        app.threads[userid] = threading.Thread(target=RUN_TREETIME, args=(root,ss))
        app.threads[userid].start()
        return jsonify({'res':'OK', 'message':'You can redirect to the wait page'})

    # error: server did not send us config for treetime run
    else:
        return jsonify({'res':'error',
            'message': 'Client-server error: server cannot find proper '
                        'config in the request'})

@app.route('/treetime/<userid>/progress', methods=['GET'])
def render_treetime_progress(userid):
    """
    Render TreeTime_progress page.
    """
    if request.method == 'GET':
        return render_template('progress_treetime.html', UserId=userid, html_theme=html_theme_css)

@app.route('/treetime/<userid>/get_session_state', methods=['GET'])
def get_session_state(userid):
    """
    In progress page, user polls the state of the session. Read the session state
    file, return the contents to the user.

    NOTE: the session file is created from another thread. Therefore, not to deal
    with asynch calls, if the file is not yet there, believe the thread will start
    normally and will cretae the file, so just return 'running' state.
    """
    root = os.path.join (sessions_root, userid)
    inf = os.path.join(root, "session_state.txt")
    if not os.path.exists(inf):
        json_data = {"session_state":{"state":"running"}}
    else:
        with open (inf, 'r') as infile:
            json_data = json.load(infile)
    return jsonify(**{"session_state": json_data})

@app.route('/treetime/<userid>/results', methods=['GET'])
def render_treetime_results(userid):
    """
    Render the treetime results page
    """
    if request.method == 'GET':
        return render_template('results_treetime.html', UserId=userid, html_theme=html_theme_css)


@app.route('/sessions/<userid>/<filename>', methods=['GET', 'POST'])
def send_file(userid, filename):
    """
    Allow user to download the requested file. The download is possible only from
    the user's session folder, i.e. 'sessions/<userid>'
    """
    uploads = os.path.join(sessions_root, userid)
    return send_from_directory(uploads, filename) #with open(os.path.join(uploads, filename), 'r') as inf:

@app.route("/ancestral/<userid>/run", methods=['POST'])
def run_ancestral(userid):
    """
    server response on the user's press 'RunAncestral' button. Reads the config
    from the html POST request, saves the config to the user's session folder.
    After then, starts a new thread and runs there the ancestral reconstruction
    job.
    """
    if (request.method != "POST"):
        abort(404)

    # save settings
    root = os.path.join(sessions_root, userid)
    if not os.path.exists(root):
        os.makedirs(root)

    #save config, run treetime in a separate thread
    if 'config' in request.get_json():
        config = request.get_json()['config']

        with open(os.path.join(root, "config.json"), 'w') as of:
            json.dump(config, of, True)
        app.threads[userid] = threading.Thread(target=RUN_TREEANC, args=(root,config))
        app.threads[userid].start()
        return jsonify({'res':'OK', 'mesage':'You can redirect to the wait page'})

    # error: server did not send us config for treetime run
    else:
        return jsonify({'res':'error',
            'message': 'Client-server error: server cannot find proper '
                        'config in the request'})

@app.route("/ancestral/<userid>/progress", methods=['GET', 'POST'])
def render_ancestral_progress(userid):
    """
    Render ancestral reconstruction progress page.
    """
    return render_template('progress_ancestral.html', UserId=userid, html_theme=html_theme_css)

@app.route("/about", methods=['GET', 'POST'])
def render_about_page():
    """
    Render about page.
    """
    return render_template('about.html', html_theme=html_theme_css)

if __name__ == "__main__":
    app.wait_time = {};
    app.threads = {};
    app.debug=True
    app.run(port=4100)
