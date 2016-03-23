import React from 'react';  
import {transitionTo, Router, browserHistory, 
  DefaultRoute, Link, Route, RouteHandler } from 'react-router';
import ReactDOM from 'react-dom'
var request = require('superagent');

import Login from './components/login.js';

var App = React.createClass({  
  render() {
    return(
      <div>
       <h1>Welcome to App</h1>
        
        <input type="button" onClick={function(){browserHistory.push('/sdrh/login');}} />
       {this.props.children}
       </div>
    );
  }
});

var Main  = React.createClass( {
  
  getInitialState() {
      
      return {
        UID: "sdf"      
      };

  },

  componentDidMount() {
    this.req_uid();
    //browserHistory.push('/zvsjdlgz/app'); 
  },

  on_uid_received(req, res, err){
    var uid = JSON.parse(res.text).redirect;
    this.setState({UID:uid});
    console.log(this.state.UID);
    browserHistory.push(this.state.UID + '/app/'); 
  },

  req_uid() {
      request
      .post('/')
      .end(this.on_uid_received);
    
  },

  render(){
    return (
        <div>
        <Link to={'/' + this.state.UID + '/app/'}>Home</Link>
        <Link to={'/' + this.state.UID + '/login/'}>Login</Link>
        {this.props.children}
        </div>
    );
  }

});

ReactDOM.render((
  <Router history={browserHistory} >
    <Route path="/" component={Main}>
      <Route path="/:user_id/app/" component={App} />
      <Route path="/:user_id/login/" component={Login} />
    </Route>
  </Router>),
document.getElementById('react'));
