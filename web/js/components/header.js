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
            <Nav>
                <NavItem eventKey={1} href="/"><Glyphicon glyph={"home"}/></NavItem>
            </Nav>
            <Navbar.Header>
              <Navbar.Brand>
                <span> TreeTime: molecular clock phylogenies</span>
              </Navbar.Brand>
              <Navbar.Toggle />
            </Navbar.Header>
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
