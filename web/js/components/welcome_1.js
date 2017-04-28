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

    on_ancestral_reconstruction_request: function(){

    },

    render:function(){
        return (
            <div>
                <Header/>
                <Button bsStyle="primary" onClick={this.ancestral_reconstruction_request}>Ancestral state reconstruction</Button>
                <Button bsStyle="primary">TreeTime run</Button>
            </div>
        );
    }
});

ReactDOM.render((
    <WelcomePage/>),
    document.getElementById('react'));

export default WelcomePage;
