from flask import Flask, Response, request,render_template,redirect, url_for, session, send_from_directory, jsonify
from werkzeug import secure_filename
import threading
import numpy as np
import os, random, subprocess, json
app = Flask(__name__)
ALLOWED_EXTENSIONS = ['fasta', 'nwk', 'csv', 'png', 'jpg']

# these modules needed to probe the input data in-place
from Bio import Phylo
import os

def make_id():
    return "HUGIYH8LJPIHP"

@app.route('/', methods=['GET', 'POST'])
def create_new_session():
    user_id = make_id();
    #return  send_from_directory('./', 'index.html')
    #return jsonify({"redirect" : user_id})
    #app.send_static_file("/templates/index.html")
    return render_template("/index.html")

@app.route('/<user_id>/', methods=['GET', 'POST'])
def open_new_session(user_id):
    return  send_from_directory('/html/index.html', './')
    #render_template("./index.html")
    




@app.route('/api/comments', methods=['GET', 'POST'])
def comments_handler():

    with open('comments', 'r') as inf:
        comments = json.loads(inf.read())

    if request.method == 'POST':
        newComment = request.form.to_dict()
        newComment['id'] = int(time.time() * 1000)
        comments.append(newComment)

        with open('comments.json', 'w') as file:
            file.write(json.dumps(comments, indent=4, separators=(',', ': ')))

    return Response(json.dumps(comments), mimetype='application/json', headers={'Cache-Control': 'no-cache', 'Access-Control-Allow-Origin': '*'})

@app.route('/probe/probe_tree', methods=['GET', 'POST'])
def probe_tree():
    if request.method == 'POST':
        print "Server probing the input newcik tree"
        data = request.form.to_dict()
        import ipdb; ipdb.set_trace()
        should_build_tree = data['should_build_tree']
        print should_build_tree
        infile = data["tree_file"]
        
        if should_build_tree=="true":
            print ("Tre will be built by the server, return OK")
            response =  jsonify({"code" : 200, "message" : "Tree will be built with FastTree method"})
            response.status_code = 200
            return response

        elif infile == "": # check file exists
            print ("Tree file is empty!")
            response =  jsonify({"code" : 500, "message" : "No tree file supplied."})
            response.status_code = 500
            return response    
        else:
        
            #try:
                tree_file = secure_filename(infile)
                tree_file.save("tree.nwk")
                with open("tree.nwk") as ff:
                    lines =  ff.readlines()
                print "lines: " + lines

                tree = Phylo.read("tree.nwk", 'newick')
                response =  jsonify({"code" : 200, "message" : "Tree loaded successfully. " + str(len(tree.get_terminals())) + " leaves."})
                response.status_code = 200
                print ("Tree file is OK")
                return response

            #except:
                print ("Tree file is corrupted.")
                response =  jsonify({"code" : 500, "message" : "Exception caught when instantiating the "
                    "tree\nMake sure you ssupopolied a valid tree in newick format."})
                response.status_code = 500
                return response

@app.route("/results/<username>")
def results(username):
    return render_template('results.html', username=username)
    #return "Hello World!"

def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1] in ALLOWED_EXTENSIONS

@app.route('/', methods=['GET', 'POST'])
def upload_file():
    if request.method == 'POST':

        print "POST REQUEST"
        data = request.form.to_dict()
        
        print "tree_file: " +  data['tree_file']

        should_build_tree = data['should_build_tree']
        tree_file = data["tree_file"] 
        aln_file = data["aln_file"] 
        meta_file = data["meta_file"]
        meta_numdate_col = data["meta_numdate_col"]
        should_build_tree = data["should_build_tree"]
        should_use_branch_len = data["should_use_branch_len"]
        should_use_branch_penalty = data["should_use_branch_penalty"] 
        branch_penalty = data["branch_penalty"]
        should_use_slope = data["should_use_slope"]
        slope_value = data["slope_value"]
        gtr = data["gtr"] 
        do_resolve_polytomies = data["do_resolve_polytomies"]
        do_coalescence_model = data["do_coalescence_model"]
        coalescence_tc = data["coalescence_tc"]
        do_autocorr_mutation_rate = data["do_autocorr_mutation_rate"]
        autocorr_alpha = data["autocorr_alpha"]
        autocorr_beta = data["autocorr_beta"]
        do_root_variance = data["do_root_variance"]

        print (should_build_tree)

        
        
        username = "".join([chr(random.randint(65,90)) for ii in range(12)])
        session['username'] = username
        user_upload_folder = os.path.join('sessions', session['username'])
        os.mkdir(user_upload_folder)
        print("all items:", request.get_json(force=True).items())
        

        print("Aligned:", request.form["aligned"])
        print("username:", username)

        seq_file = request.files['sequence_file']
        if seq_file and allowed_file(seq_file.filename):
            filename = secure_filename(seq_file.filename)
            seq_file.save(os.path.join(user_upload_folder, "sequences.fasta"))
            #print url_for('uploaded_file', filename=filename)
        tree_file = request.files['tree_file']
        if tree_file and allowed_file(tree_file.filename):
            filename = secure_filename(tree_file.filename)
            tree_file.save(os.path.join(user_upload_folder, "initial_tree.nwk"))
            #print url_for('uploaded_file', filename=filename)
        meta_data_file = request.files['meta_data_file']
        if meta_data_file and allowed_file(meta_data_file.filename):
            filename = secure_filename(meta_data_file.filename)
            meta_data_file.save(os.path.join(user_upload_folder, "meta_data.csv"))
            #print url_for('uploaded_file', filename=filename)

        app.threads[username] = threading.Thread(target = process, args = (username,))
        app.threads[username].start()
        app.wait_time[username] = 1
        return redirect("progress/"+username)

    #return render_template('welcome.html')
    return render_template("welcome3.html")

@app.route("/progress/<username>")
def progress(username):
    app.wait_time[username] = min(2*app.wait_time[username], 10)
    steps = [
                ["todo", "building tree using fasttree"],
                ["todo", "branch length optimzation"],
                ["todo", "date constraints"],
                ["todo", "estimate internal node locations"],
                ["todo", "export"]]
    try:
        with open(os.path.join('sessions', username, 'status.log')) as logfile:
            log = logfile.readlines()
        steps_done = [s.strip() for s in log]
        for step in steps:
            if step[1] in steps_done:
                step[0]="done"
    except:
        log = ["status file opening failed"]

    if log[-1].strip()=='done':
        print("successfully ran the analysis, redirecting ",username, " to results")
        return redirect('/results/'+username)
    elif app.wait_time[username]>1000:
        print("TIMEOUT, redirecting ",username, " to results")
        return redirect('/results/'+username)
    else:
        print(steps)
        return render_template("wait.html", wait=app.wait_time[username], steps=steps)

@app.route("/static/")
def treetime_static():
    return url_for('static', filename='css/style.css')

@app.route("/sessions/<path:path>")
def send_user(path):
    return send_from_directory('sessions', path)



if __name__ == "__main__":
    app.wait_time = {};
    app.threads = {};
    tmp_key = os.urandom(24)
    app.secret_key = tmp_key
    app.debug=True
    app.run(port=3000)
