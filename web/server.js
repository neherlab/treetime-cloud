var express = require('express');
var app = express();
var multer = require("multer");
var fs = require('fs');
var bodyParser = require('body-parser')
var upload = multer({ dest: 'uploads/' })

var subp = require('child_process')

app.use(express.static(__dirname ));
// parse application/json
app.use(bodyParser.json())

function makeid()
{
    var text = "";
    var possible = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";

    for( var i=0; i < 7; i++ )
        text += possible.charAt(Math.floor(Math.random() * possible.length));

    return text;
}

var mkdirSync = function (path) {
  try {
    fs.mkdirSync(path);
    console.log("Session directory created!")
  } catch(e) {
    if ( e.code != 'EEXIST' ) throw e;
  }
}

app.get('/', function (req, res) {
  var user_id = makeid(); //"HILWQ89P23OI566UHLUIL";
  res.send(user_id);
});

app.post('/', function(req, res){
  
  var user_id = makeid(); //var user_id = "HILWQ89P23OI566UHLUIL";
  console.log('User requests for a new ID, send' + user_id + " back")
  mkdirSync(__dirname + "/sessions/" + user_id);  
  res.send({redirect:user_id});
  

});


app.get('/:user_id/', function (req, res) {
  
  res.sendfile( __dirname + "/static/" + 'index.html');

});


app.post('/:user_id/tree_file', upload.single('treefile'), function(req, res){
  var user_id = req.params.user_id;
  console.log("USER_ID: " + user_id);
  fs.rename(__dirname + "/" + req.file.path, __dirname + "/sessions/" + user_id + "/in_tree.nwk", function (err) {
  console.log('ERR: ');
  console.log(err);
});
  console.log(req.file)
});

app.post('/:user_id/aln_file', upload.single('alnfile'), function(req, res){
  var user_id = req.params.user_id;
  console.log("USER_ID: " + user_id);
  fs.rename(__dirname + "/" + req.file.path, __dirname + "/sessions/" + user_id + "/in_aln.fasta", function (err) {
  console.log('ERR: ');
  console.log(err);
});
  console.log(req.file)
});

app.post('/:user_id/meta_file', upload.single('metafile'), function(req, res){
  var user_id = req.params.user_id;
  console.log("USER_ID: " + user_id);
  fs.rename(__dirname + "/" + req.file.path, __dirname + "/sessions/" + user_id + "/in_meta.csv", function (err) {
  console.log('ERR: ');
  console.log(err);
});
  console.log(req.file)
});


app.post('/:user_id/run',  function(req, res){

  console.log(req.params.user_id);
  var user_id = req.params.user_id


  fs.writeFile(__dirname + "/sessions/" + user_id + "/settings.json", JSON.stringify(req.body.settings), function(err) {
    if(err) {
      res.send({status:"Error", message:"Cannot save ssession settings"});
    }
  }); 

  var root = __dirname + "/sessions/" + user_id;
  subp.exec_sync('python ./static/py/main.py ' + root, function(err, stdout, stderr){

  });


  res.send({status:"OK"})

});

app.post('/:user_id/session_state', function(req, res){
  var user_id = req.params.user_id
  console.log(user_id)
  var fname = __dirname + "/sessions/" + user_id + "/state.json"
  fs.readFile(fname, "utf8", function(err, data) {
        if (err) throw err;
        res.send(data);
    });
});
 
app.post('/:user_id/results', function(req, res){
  var user_id = req.params.user_id
  console.log(user_id)
  var fname = __dirname + "/sessions/" + user_id + "/tree.json"
  fs.readFile(fname, "utf8", function(err, data) {
        if (err) throw err;
        res.send({tree:JSON.parse(data)});
    });
});
 

app.listen(3000, function () {
  console.log('Example app listening on port 3000!');
});
