import React from  'react'
import ReactDOM from 'react-dom'
import Header from './header.js'
var request = require('superagent');
import { Button } from "react-bootstrap";

var WelcomeAncPage = React.createClass({

    componentDidMount: function(){
        var parentNode = this.getDOMNode().parentNode;
        // set user id from the HTML template
        this.setAppState({'UID':parentNode.attributes.userid.value});
    },

    setAppState : function(partialState){
        this.setState(partialState);
    },

    render:function(){
        return (
            <div>
                <Header/>
                <Button bsStyle="primary" onClick={}>Alignment</Button>
                <Button bsStyle="primary">Tree</Button>
            </div>
        );
    }
});

ReactDOM.render((
    <WelcomeAncPage/>),
    document.getElementById('react'));

export default WelcomeAncPage;
