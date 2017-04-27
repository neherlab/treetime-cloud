import React from  'react'
import ReactDOM from 'react-dom'
import Header from './header.js'

var WelcomePage = React.createClass({
    render:function(){
        return (

            <div>
                <Header/>
            </div>
        );
    }
});

ReactDOM.render((
    <WelcomePage/>),
    document.getElementById('react'));

export default WelcomePage;
