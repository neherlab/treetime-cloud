import React from  'react'
import ReactDOM from 'react-dom'
import Header from './header.js'
var request = require('superagent');
import { Button } from "react-bootstrap";

var WelcomePage = React.createClass({

    componentDidMount: function(){

    },

    setAppState : function(partialState){
        this.setState(partialState);
    },

    ancestral_reconstruction_request: function(){

        request.post("/ancestral_reconstruction_request")
              .set('Content-Type', 'application/json')
              //.send({config: this.state.config})
              .end(this.on_ancestral_reconstruction_request);
    },

    on_ancestral_reconstruction_request: function(err, res){

        if (err){
            alert("Cannot perform ancestral state reconstruction: Internal server error.")
            return;
        }
        if (res.statusCode != 200){
            alert("Cannot perform ancestral state reconstruction: Internal server error.\nServer says: " + res.statusText)
            return;
        }
        if (!res.text){
            alert("Cannot perform ancestral state reconstruction: Server unable to start new session.")
            return;
        }
        var data = JSON.parse(res.text)
        var uid = data.UserId;
        if (!uid){
            alert("Cannot perform ancestral state reconstruction: Server unable to assign new user id.")
            return;
        }
        this.setAppState({'UID':uid});
        window.location.replace("/ancestral/" + data.UserId);

    },

    treetime_request: function(){

        request.post("/treetime_request")
              .set('Content-Type', 'application/json')
              //.send({config: this.state.config})
              .end(this.on_treetime_request);
    },

    on_treetime_request: function(err, res){

        if (err){
            alert("Cannot start treetime session: Internal server error.")
            return;
        }
        if (res.statusCode != 200){
            alert("Cannot start new treetime session: Internal server error.\nServer says: " + res.statusText)
            return;
        }
        if (!res.text){
            alert("Cannot start new treetime session: Server did not provide enough data.")
            return;
        }
        var data = JSON.parse(res.text)
        var uid = data.UserId;
        if (!uid){
            alert("Cannot start new treetime session: Server unable to assign new user id.")
            return;
        }
        this.setAppState({'UID':uid});
        window.location.replace("/treetime/" + data.UserId);
    },

    render:function(){
        return (
            <div>
                <Header/>
                <Button bsStyle="primary" onClick={this.ancestral_reconstruction_request}>Ancestral state reconstruction</Button>
                <Button bsStyle="primary" onClick={this.treetime_request}>TreeTime run</Button>
            </div>
        );
    }
});

ReactDOM.render((
    <WelcomePage/>),
    document.getElementById('react'));

export default WelcomePage;
