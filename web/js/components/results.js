import React from  'react'
var request = require('superagent');

var Res = React.createClass({
    
    getInitialState(){
        return (
        {
            tree:{} 
        }
        );
    },

    componentWillMount(){
        request.post('/' + this.props.UID + '/results')
            .send({UID:this.state.UID})
            .end(this.on_results);
    },

    on_results(err, res){
        if (err) return;

        var data = JSON.parse(res.text);
        if (data.tree){
            this.setState({tree:data.tree});    
        }
    },

    render(){

        return (
            <div>
            <h1>Results</h1>
            </div>
        );
    }
});

export default Res;