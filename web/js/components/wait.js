import React from  'react'
var request = require('superagent');

var Banner = React.createClass({
    render(){
        return (
            <div  id="waitbanner" >
            <div style={{"width":"100%"}}>
                <h3>TreeTime is running...</h3>
                <h4>This may take a while. You will be redirected to the results page when its done</h4>
            </div>
            </div>
            );
    }
});

var Status = React.createClass({
    render(){
        return (
            <div>
                <h4>{this.props.status}</h4>
            </div>
        );
    }
});

var StepName = React.createClass({
    render(){
        return (

            <div>
            
            <h4>{this.props.step}</h4>
            </div>
        );
    }
});

var Step = React.createClass({
    
    render() {
        return(
            <div>
                <Status status={this.props.status}/>
                <StepName step={this.props.step}/>
            </div>
        );
    }
});

var Wait  = React.createClass({
    getInitialState(){
        return (
        {
            UID : "Undef",
            steps: []
        }
        );
    },

    componentWillMount(){
        console.log(this.props.UID)
        this.state.user_id = this.props.UID
        console.log("Main: user_id = " + this.state.user_id)
        //setInterval(this.get_tree_res, 2000);
    },

    request_session_state(){
       
        request.post('/' + this.props.UID + '/session_state')
            .send({UID:this.state.UID})
            .end(this.on_session_state);
    },

    on_session_state(err, res){
        console.log(JSON.parse(res.text));
        //var sts = JSON.parse(res);
        this.setState({steps:JSON.parse(res.text)});
        console.log(this.state.steps)
        if (this.check_all_done()){
            clearInterval(this.interval);
            this.props.handle_all_done();
        }
    },

// /    post_callback:(err, res, body) ->
// /        
// /        console.log(res);
// /        if err
// /            console.log ("Run request failed: " + err);
// /        else
// /            data = JSON.parse(body);
// /            this.setState({steps: data.steps});
// /            
// /            arr = Object.keys(this.state.steps).map(this.check_done)
// /            if not false in arr
// /                console.log("ALL done!")
// /                clearInterval(this.interval);
// /
    check_all_done(){
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



    componentDidMount(){
        this.request_session_state();
        this.interval = setInterval(this.request_session_state, 2000);
    },
    
    componentWillUnmount(){
        clearInterval(this.interval);
    },
    
    render_step (step){
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

    render(){
        
        return (
            <div>
                <Banner/>
                {(this.state.steps).map(this.render_step)}
            </div>
        );
    }
});


export default Wait;

