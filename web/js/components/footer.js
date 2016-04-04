import React from 'react'

var Copyright = React.createClass({
    render(){
        var scope = {
            splitterStyle: {
                width: 100
            }
        };
        return (
            <div id='copyright'>
                Copyright &copy;:
            </div>

        );
    }
});

var Authors = React.createClass({

    render(){
        return (
            <div id="authors">
            Richard Neher and Pavel Sagulenko, 2016
            </div>

        );
    }
});

var Footer  = React.createClass({

    render (){
        return (
            <div id="footer">
                <Authors/>
                <Copyright/>
            </div>
        );
    }

});

export default Footer; 
