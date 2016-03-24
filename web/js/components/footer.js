import React from 'react'

var Copyright = React.createClass({
    render(){
        var scope = {
            splitterStyle: {
                width: 100
            }
        };
        return (
            <div>
                Copyright:
            </div>

        );
    }
});

var Authors = React.createClass({

    render(){
        return (
            <div>
            Richard Neher and Pavel Sagulenko, 2016
            </div>

        );
    }
});

var Footer  = React.createClass({

    render (){
        return (
            <div id="footer">
                <Copyright/>
                <Authors/>
            </div>
        );
    }

});

export default Footer; 
