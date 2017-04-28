import React from  'react'
import ReactDOM from 'react-dom'
import Header from './header.js'
var request = require('superagent');
import { Button } from "react-bootstrap";

var WelcomeTreeTimePage = React.createClass({

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
                <Button bsStyle="primary">TreeTime</Button>
                <Button bsStyle="primary">TreeTime</Button>
            </div>
        );
    }
});

ReactDOM.render((
    <WelcomeTreeTimePage/>),
    document.getElementById('react'));

export default WelcomeTreeTimePage;
