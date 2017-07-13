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
            Pavel Sagulenko and Richard Neher, <a href="https://doi.org/10.1101/153494">bioRxiv</a>, 2017.
            </div>
        );
    }
});

var Footer  = React.createClass({

    render (){
        return (
            <div id="footer_wrapper" style={{"position":"absolute",  "bottom":"0%", "left":"0"}}>
                <a href='/about' id="impressum">About/Impressum</a>
                <div id="footer">
                    <Authors/>
                    <Copyright/>
                </div>
            </div>
        );
    }

});

export default Footer;
