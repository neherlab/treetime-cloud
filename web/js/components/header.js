import React from 'react'


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
            <h1>TimeTree </h1>
            <h2>Phylogeny with absolute time resolution</h2>
            </div>

        );
    }
});

var Header  = React.createClass({

    render(){
        return (
            <div id="header">
                <Logo/>
                <Name/>
            </div>
        );
    }
});

export default Header;