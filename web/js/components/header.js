import React from 'react'
import { Row } from "react-bootstrap";


var Logo = React.createClass({

    render(){
        var scope = {
            splitterStyle: {
                width: 100
            }
        };
        return (
                <img id='logo' src="/static/svg/logo.svg" />
        );
    }

});

var Name = React.createClass({

    render(){
        return (
            <div id='name'>
            <h1>TreeTime </h1>
            <div className="spacer"></div>
            <h4>Maximum likelihood Inference molecular clock phylogenies</h4>
            </div>
        );
    }
});

var Header  = React.createClass({

    render(){
        return (
            <Row>
            <div id="header">
                <Logo/>
                <Name/>
            </div>
            </Row>
        );
    }
});

export default Header;