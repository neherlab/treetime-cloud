import React from 'react'


var Logo = React.createClass({

    render(){
        var scope = {
            splitterStyle: {
                width: 100
            }
        };
        return (
            <div id='logo'>
                <h1>Logo</h1>
            </div>

        );
    }

});

var Name = React.createClass({

    render(){
        return (
            <div id='name'>
            <h1>TimeTree </h1>
            <h2>Phylogeny with dating</h2>
            </div>

        );
    }
});

var Header  = React.createClass({

    render(){
        return (
            <div id="header" style={{"width":"100%"}}>
                <Logo/>
                <Name/>
            </div>
        );
    }
});

export default Header;