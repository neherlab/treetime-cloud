import React from 'react';  
import {transitionTo, Router, browserHistory, 
  DefaultRoute, Link, Route, RouteHandler } from 'react-router';
import ReactDOM from 'react-dom'
var request = require('superagent');

import Login from './components/login.js';
import Header from './components/header.js'
import Footer from './components/footer.js'
import TreeTimeForm from './components/welcome.js'
import Wait from './components/wait.js'
import Res from "./components/results.js"
var settings = {
  
  doBuildTree:true,
  shouldReuseBranchLen:false,
  doReroot:false,
  gtr:"Jukes-Cantor",
  shouldUseBranchLenPenalty:{
      bool:true,
      value:0.0
  },
  shouldUseSlope:{
      bool:true,
      value:0.0
  },
  doResolvePoly: false,
  doCoalescent:{
    bool:true,
    Tc:0.0
  },
  doRelaxedClock:{
    bool:true,
    alpha:0.0,
    beta:0.0
  },
  doRootJoint:true,
  doCalcGTR:true
};

var Main  = React.createClass( {
  
  getInitialState() {
      return {
        UID: "sdf", 
        settings:settings,
        state: {},      
      };
  },
  
  componentDidMount() {
    this.req_uid();
    //browserHistory.push('/zvsjdlgz/app'); 
  },
  
  req_uid() {
      request
      .post('/')
      .end(this.on_uid_received);
    
  },
  
  on_uid_received(req, res, err){
    var uid = JSON.parse(res.text).redirect;
    this.setState({UID:uid});
    console.log(this.state.UID);
    browserHistory.push(this.state.UID + '/app/'); 
  },
  
  handle_run(){
    console.log("APP:: RUN button pressed");
    request.post("/" + this.state.UID + "/run/")
      .set('Content-Type', 'application/json')
      .send({settings: this.state.settings})
      .end(this.on_run_status);
  },

  on_run_status(err, res){
    var status = JSON.parse(res.text).status;
    console.log(status);
    if (status == "OK"){
      browserHistory.push('wait/'); 
    }
  },

  on_settings_changed(name, setting){
      console.log("APP:: settings changed. " + name + " new value = " + setting);
      var settings = this.state.settings
      settings[name] = setting;

      this.setState({settings: settings})
      //this.state.settings = settings;
  },

  on_state_changed(name, state){

  },

  on_all_done(){
      console.log("ALL don, redirecting to RESULTS page");
      browserHistory.push('results/'); 
  },

  render(){
    return (
        <div>
        <Header/>
        {this.props.children && 
          React.cloneElement(
            this.props.children, 
            {
              UID:this.state.UID,
              settings:this.state.settings,
              state:this.state.state,
              handle_run: this.handle_run,
              handle_settings_change: this.on_settings_changed,
              handle_state_changed: this.on_state_changed, 
              handle_all_done: this.on_all_done
            }
          )
        }
        <Footer/>
        </div>
    );
  }
});

ReactDOM.render((
  <Router history={browserHistory} >
    <Route path="/" component={Main}>
      <Route path="/:user_id/app/" component={TreeTimeForm} />
      <Route path="wait/" component={Wait} />
      <Route path="results/" component={Res} />
    </Route>
  </Router>),
document.getElementById('react'));
