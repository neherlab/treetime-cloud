var express = require('express');
var app = express();

app.use(express.static(__dirname ));

app.get('/', function (req, res) {
  var user_id = "HILWQ89P23OI566UHLUIL";
  res.send(user_id);
  //res.redirect('http://google.com');
  //res.redirect('./' + user_id + "/");
  //res.sendfile( 'index.html');
});

app.post('/', function(req, res){
  
  var user_id = "HILWQ89P23OI566UHLUIL";
  console.log('User requests for a new ID, send' + user_id + " back")
  res.send({redirect:user_id});
  

});

app.get('/:user_id/', function (req, res) {
  res.sendfile( 'index.html');
});

app.listen(3000, function () {
  console.log('Example app listening on port 3000!');
});
