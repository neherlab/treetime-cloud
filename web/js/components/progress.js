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

var Status = React.createClass({

    render(){
        return (
            <div>
                <h4>{this.props.status}</h4>
            </div>
        );
    }

});

var StepName = React.createClass({

    render(){
        return (

            <div>
            <h4>{this.props.step}</h4>
            </div>
        );
    }
});

var Step = React.createClass({

    style:{
        'float': 'left',
        'width':'100%',
        'posiiton':'relative',
        'text-align':'left',
    },

    render() {
        return(
            <div style={this.style}>
                <Status status={this.props.status}/>
                <StepName step={this.props.step}/>
            </div>
        );
    }
});

var Progress  = React.createClass({

    getInitialState(){
        return (
        {
            UID : "Undef",
            steps:[],
            todos : [],
            progress : "",
            done : [],
            session_state : "None",
            error : "",
            warns : []
        }
        );
    },

    request_session_state: function(){

        request.get('/' + this.state.UID + '/session_state')
            .send({UID:this.state.UID})
            .end(this.on_session_state);
    },

    on_session_state : function (err, res){
        if (err){
            if (err.status == 404){
                window.location.replace("404")
            }
            return;
        }

        var state = JSON.parse(res.text).steps
        this.setState({
            done:state.done,
            error: state.error,
            progress: state.progress,
            session_state:state.state,
            todos:state.todo,
            warns:state.warn
        });

        if (this.state.error != ""){
            clearInterval(this.interval);
        }
        console.log("Check session state...")
        if (this.check_all_done()){
            console.log("All Done!")
            clearInterval(this.interval);
            window.location.replace('/' + this.state.UID + '/results')
        }
    },

    componentDidMount: function(){

        var parentNode = this.getDOMNode().parentNode;
        //console.log(parentNode)
        var UID = (parentNode.attributes.userid.value);
        //console.log("UID: " + UID)
        this.state.UID = UID;
        this.request_session_state();
        this.interval = setInterval(this.request_session_state, 10000);
    },

    setAppState : function(partialState){
        this.setState(partialState);
    },

    componentWillUnmount : function(){
        clearInterval(this.interval);
    },


    check_all_done : function(){
        console.log(this.state.session_state);
        if (this.state.session_state == "complete"){
            return true;
        }else{
            return false;
        }
    },


    render: function(){

        return (
            <div>
                <Header />
                <Body
                    appState={this.state}
                    setAppState={this.setAppState}/>
                <Footer />
            </div>
        );
    }
});

var Body = React.createClass({

    render_step : function (type, step){

        // pick a glyphicon from bootstrap
        var glyph_type = "";
        switch(type){

            case ("todo"):
            glyph_type = "option-horizontal";
            break;

            case ("progress"):
            glyph_type = "arrow-right";
            break;

            case ("done"):
            glyph_type = "ok";
            break;
        }



        //console.log("status = "+ s)
        if (status != 'Error'){
            return (
            <Row className="grid-treetime-row">
                <Col xs={1} md={1} id="results_col_step" className="grid-treetime-col-left">
                    <Glyphicon glyph={glyph_type}/>
                </Col>

                <Col xs={11} md={11} id="results_col_step" className="grid-treetime-col-right">
                    <span>{step}</span>
                </Col>

            </Row>
            );
        }
        else {
            // TODO show error banner
            clearInterval(this.interval);
        }
    },

    render_error: function(error){
        var glyph_type = "ban-circle"
        var mail = "mailto:pavel.sagulenko@tuebingen.mpg.de?Subject=TreeTime%20error. Session:" + this.props.appState.UID
        return (
            <div className="progress_error">
            <h4 id="error_header"> Ooops... An error occured</h4>
            <Panel collapsible defaultExpanded header="" id='error_message'>
            <Row className="grid-treetime-row">
                <Glyphicon id="error_header" glyph={glyph_type}/>
            </Row>
            <Row className="grid-treetime-row">
                <span id="error_header" style={{"font-weight": "bold"}}>Server says: {error}</span>
            </Row>
            <Row className="grid-treetime-row">
                <p style={{"text-align":"justify"}}>
                Either the input parameters caused some numerical problems/overflow exceptions, or you encountered a bug in TreeTime.
                Please, send us an <a href={mail} target="_top">e-mail</a> to help us diagnose and fix the problem. (please keep the session ID in the mail subject). THANK YOU.
                </p>
            </Row>

            <Row className="grid-treetime-row">
                <p style={{"text-align":"justify"}}>
                We are sorry for the incovenience and will try to fix the problem as soon as possible.
                </p>
            </Row>

            </Panel>
            </div>
        )
    },

    switch_error: function(){
        if (this.props.appState.error != ""){
            return this.render_error(this.props.appState.error)
        }else{
            return null;
        }
    },

    render: function(){
        var render_step = this.render_step
        return (
            <div className="page_container">
                <Banner/>
                <div className="bigspacer"></div>
                <Panel header="Progress...">
                {(this.props.appState.done).map(function(d){
                    return render_step("done", d)
                })}

                {render_step("progress", this.props.appState.progress)}

                {this.switch_error()}

                {(this.props.appState.todos).map(function(d){
                    return render_step("todo", d)
                })}
                </Panel>

            </div>
            );


    }
})

ReactDOM.render((
    <Progress />),
    document.getElementById('react'));

export default Progress;

