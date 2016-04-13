import React from  'react'
import ReactDOM from 'react-dom'
import Header from './header.js'
import Footer from './footer.js'

var request = require('superagent');

var Banner = React.createClass({
    render(){
        return (
            <div>
            <div>
                <h3>TreeTime is running...</h3>
                <h4>This may take a while. You will be redirected to the results page when its done</h4>
            </div>
            </div>
            );
    }
});

var Status = React.createClass({

    statusStyle : {
        'float': 'left',
        'width':'10%',
        //'border-style' : 'solid',
        //'border-width' : '1px'
    },

    render(){
        return (
            <div style={this.statusStyle}>
                <h4>{this.props.status}</h4>
            </div>
        );
    }
});

var StepName = React.createClass({
    
    style : {
        'float': 'left',
        'width':'80%',
        //'border-style' : 'solid',
        //'border-width' : '1px'
    },

    render(){
        return (

            <div style={this.style}>
            <h4>{this.props.step}</h4>
            </div>
        );
    }
});

var Step = React.createClass({
    
    style:{
        'float': 'left',
        'width':'100%',
        'posiiton':'relative',
        'text-align':'left',
    },

    render() {
        return(
            <div style={this.style}>
                <Status status={this.props.status}/>
                <StepName step={this.props.step}/>
            </div>
        );
    }
});

var Progress  = React.createClass({
    
    getInitialState(){
        return (
        {
            UID : "Undef",
            steps: []
        }
        );
    },

    request_session_state: function(){
       
        request.get('/' + this.state.UID + '/session_state')
            .send({UID:this.state.UID})
            .end(this.on_session_state);
    },

    on_session_state : function (err, res){
        if (err){
            console.log(err);
            console.log(err.status)
            if (err.status == 404){
                window.location.replace("404")
            }
            return;
        }
        console.log(res.body)
        //console.log(JSON.parse(res.text));
        ////var sts = JSON.parse(res);
        this.setState({steps:JSON.parse(res.text).steps});
        console.log("Steps:")
        console.log(this.state.steps)
        console.log(JSON.parse(res.text))
        if (this.check_all_done()){
            clearInterval(this.interval);
            window.location.replace('/' + this.state.UID + '/results')
            // TODO redirect
        }
    },

    componentDidMount: function(){
        
        var parentNode = this.getDOMNode().parentNode;
        console.log(parentNode)
        var UID = (parentNode.attributes.userid.value);
        console.log("UID: " + UID)
        this.state.UID = UID;
        this.request_session_state();
        this.interval = setInterval(this.request_session_state, 10000);

    },
 
    setAppState : function(partialState){
        this.setState(partialState);
    },
    
    componentWillUnmount : function(){
        clearInterval(this.interval);
    },


    check_all_done : function(){
        var arrayLength = this.state.steps.length;
        for (var i = 0; i < arrayLength; i++) {
            var step =  this.state.steps[i];
            console.log(step);
            if (step.status != "Done"){
                return false;
            }
        }
        return true;
    },

    progressStyle : {
        'float':'left',
        'position':'relative',
        'text-align':'center',
        'width':'100%',
        //'border-style':'solid',
        //'border-width': '1px',
    },

    render: function(){
        
        return (
            <div style={this.progressStyle}>
                <Header />
                <Body 
                    appState={this.state}
                    setAppState={this.setAppState}/>
                <Footer />
            </div>
        );
    }
});

var Body = React.createClass({
   
   bodyStyle : {
    'width':'70%',
    //'float':'left',
    'display': 'inline-block',
    'position':'relative',
    //'border-style':'solid',
    //'border-width': '1px'
    }, 
    
    render_step : function (step){
        console.log("rendering" + step.name)
        console.log(step)
        
        var s = step.status
        var n = step.name
        var key = step.name
        console.log("status = "+ s)
        if (status != 'Error'){
            return <Step key={key} status={s} step={n} />
        }
        else {
            // TODO show error banner 
            clearInterval(this.interval);
        }
    },

    render: function(){
        return (
            <div style={this.bodyStyle}>
                <Banner/>
                {(this.props.appState.steps).map(this.render_step)}
            </div>
        );
    }
})

ReactDOM.render((
    <Progress />),
    document.getElementById('react'));

export default Progress;

