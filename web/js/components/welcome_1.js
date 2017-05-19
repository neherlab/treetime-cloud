import React from  'react'
import ReactDOM from 'react-dom'
import Header from './header.js'
import Footer from './footer.js'
var request = require('superagent');
import { Button, Row, Col, Panel } from "react-bootstrap";

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
        var btn_size={"width":"200px", "height":"80px", "margin-right":"25px", "margin-left":"25px", "display":"inline-block", "vertical-alignment":"middle"};
        return (


            <div>
                <Header/>
                <div className="page_container">
                    <div id="welcome_description">
                        <Panel collapsible defaultExpanded header="Time-tree inference">
                        <Row>
                            <Col md={3}>
                                <div style={{"height":"30px"}}></div>
                                <Button bsStyle="primary" style={btn_size} onClick={this.treetime_request}>Timetree inference</Button>
                            </Col>

                            <Col md={5}>
                                <h4>Features</h4>
                                <ul>
                                    <li>Approximate maximum-likelihood time tree inference</li>
                                    <li>Inference of GTR models</li>
                                    <li>Rerooting to obtain best root-to-tip regression</li>
                                    <li>Coalesent priors</li>
                                    <li>Auto-correlated molecular clocks</li>
                                </ul>
                            </Col>
                            <Col md={3}>
                                <h4>Requires</h4>
                                <ul>
                                    <li>Alignment as fasta</li>
                                    <li>Tip dates as csv</li>
                                    <li>Tree as newick (optional)</li>
                                </ul>
                            </Col>
                        </Row>
                        </Panel>
                        <Panel collapsible defaultExpanded header="Ancestral state reconstruction">
                        <Row>
                            <Col md={3}>
                                <div style={{"height":"30px"}}></div>
                                <Button bsStyle="primary" style={btn_size} onClick={this.ancestral_reconstruction_request}>Ancestral sequence<br/>&#13;&#10;reconstruction</Button>
                            </Col>

                            <Col md={5}>
                                <h4>Features</h4>
                                <ul>
                                    <li>Ancestral sequence reconstruction</li>
                                    <li>Branch length optimization</li>
                                    <li>Inference of GTR models</li>
                                </ul>
                            </Col>
                            <Col md={3}>
                                <h4>Requires</h4>
                                <ul>
                                    <li>Alignment as fasta</li>
                                    <li>Tree as newick (optional)</li>
                                </ul>
                            </Col>
                        </Row>
                        </Panel>
                    </div>
                </div>
                <Footer/>
            </div>
        );
    }
});

ReactDOM.render((
    <WelcomePage/>),
    document.getElementById('react'));

export default WelcomePage;
