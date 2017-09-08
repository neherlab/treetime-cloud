import React from 'react'
import { Row, Navbar, Nav, NavItem, Glyphicon } from "react-bootstrap";

<h4>Maximum-Likelihood molecular clock phylogenies</h4>
var Logo = React.createClass({

    render(){
        var scope = {
            splitterStyle: {
                width: 30
            }
        };
        return (

                <img id='logo' style={{"width":"30px;"}} src="/static/svg/logo.svg" />
        );
    }

});

var Name = React.createClass({

    render(){
        return (
            <div id='name'>
            <h3>TreeTime </h3>
            </div>
        );
    }
});

var Header  = React.createClass({
    render(){
        const navbarInstance = (
          <Navbar>
            <Navbar.Header>
              <Navbar.Brand style={{"fontSize":"20pt"}}>
                TreeTime: Maximum-likelihood phylodynamic analysis
              </Navbar.Brand>
              <Navbar.Toggle />
            </Navbar.Header>
            <Nav pullRight>
                <NavItem eventKey={1} href="/"><Glyphicon glyph={"home"}/> Home</NavItem>
                <NavItem eventKey={2} href="/about"><Glyphicon glyph={"info-sign"}/> About</NavItem>
                <NavItem eventKey={3} href="/doc"><Glyphicon glyph={"list-alt"}/> API doc</NavItem>
            </Nav>

          </Navbar>
        );

        return (
            navbarInstance
            // <Row>
            // <div className="page_container">
            // <div id="header">
            //     <Logo/>
            //     <Name/>
            // </div>
            // </div>
            // </Row>
        );
    }
});

export default Header;
