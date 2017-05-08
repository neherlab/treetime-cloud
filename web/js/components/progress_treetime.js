import React from  'react'
import ReactDOM from 'react-dom'
import Header from './header.js'
import Footer from './footer.js'

import { Glyphicon, Panel, Button, Grid, Row, Col, FormControl} from "react-bootstrap";

var request = require('superagent');

var Banner = React.createClass({
    render(){
        return (
            <div>
            <div>
                <h3>TreeTime is running...</h3>
                <h4>This may take a while. You will be redirected to the results page when its done</h4>
            </div>
            </div>
            );
    }
});

var ErrorBanner = React.createClass({


    render(){
        var glyph_type = "ban-circle"
        var mail = "mailto:pavel.sagulenko@tuebingen.mpg.de?subject=TreeTime%20error.%20Session:" + this.props.appState.UID

        return (
            <div>
                <h3 id="error_header">Oooops... An error occured.</h3>
                <div>
                <p style={{"text-align":"justify"}}>
                    Either the input parameters caused some numerical problems/overflow exceptions, or you encountered a bug in TreeTime.
                    Please, send us an <a href={mail} target="_top">e-mail</a> to help us diagnose and fix the problem. (please keep the session ID in the mail subject). THANK YOU.
                </p>
                </div>
                <div>
                    <p style={{"text-align":"justify"}}>
                    We are sorry for the incovenience and will try to fix the problem as soon as possible.
                    </p>
                </div>
                <Panel collapsible defaultCollapsed header="Server output:" id='error_message'>
                    {this.props.appState.error}
                </Panel>
            </div>
        );
    }
})

var ProgressTreeTimePage = React.createClass({

    componentDidMount: function(){
        var parentNode = this.getDOMNode().parentNode;
        //console.log(parentNode)
        var UID = (parentNode.attributes.userid.value);
        //console.log("UID: " + UID)
        this.state.UID = UID;
        this.requestSessionSate();
        this.interval = setInterval(this.requestSessionSate, 10000);
    },

    getInitialState(){
        return (
            {
                UID : "",
                "state":true,
                "current_state":"running",
                "error":""
            }
        );
    },

    setAppState : function(partialState){
        this.setState(partialState);
        console.log("Progress page setting app state: ")
        console.log(this.state)
    },

    requestSessionSate: function(){
        console.log("Request session state sent...")
        request.get('/treetime/' + this.state.UID + '/get_session_state')
            .send({UID:this.state.UID})
            .end(this.onSessionState);
    },

    onSessionState : function (err, res){


        if (err){
            if (err.status == 404){
                window.location.replace("404")
            }
            return;
        }

        var session_state = JSON.parse(res.text).session_state
        if (!session_state){
            this.setAppState({"state":false, "error":"Internal server error: The server failed to send treetime status update."})
            return;
        }
        console.log("Current state: " + session_state.state)
        this.setAppState({"error":session_state.desc, "state":session_state.state!="error", "current_state" : session_state.state})
        if (this.checkDone()){
            window.location.replace("/treetime/" + this.state.UID + "/results");
        }
    },

    checkDone: function(){
        return this.state.current_state == "done"
    },

    myXOR : function(a,b) {
        return ( a || b ) && !( a && b );
    },

    stateRenderStyle: function(error_page){
        return this.myXOR(error_page, this.state.state) ? {"display":"inline-block"} : {"display":"none"}
    },

    render: function(){
        return (
            <div>
                <div style={this.stateRenderStyle(false)}>
                    <Banner  />
                </div>
                <div style={this.stateRenderStyle(true)}>
                    <ErrorBanner appState={this.state}/>
                </div>
            </div>
        )
    }
})

ReactDOM.render((
    <ProgressTreeTimePage />),
    document.getElementById('react'));

export default ProgressTreeTimePage;

