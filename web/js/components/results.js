import React from  'react'
var request = require('superagent');
import PhyloTree from './phylo_tree.js'
import TreeLegend from './tree_legend.js'
import MuPlot from './mu_plot.js'
import RootLhPlot from './root_lh.js'

var Globals = require('./globals.js')
var colors = Globals.colors;

var CScale = function(cUnit){

    this.colors = ["#4D92BF", "#5AA5A8", "#6BB18D", "#80B974", "#98BD5E", "#B1BD4E",
          "#C8B944", "#DAAC3D", "#E59738", "#E67732", "#E14F2A", "#DB2522"];
    
    this.cUnit = cUnit;
    
    this.get_color = function(x){
        if (!this.color){
            return colors[5];
        }else{
            return this.color(x);
        }
    };
    
    this.update = function (cUnit, root){
        if (!root) return;
        this.cUnit = cUnit;
        var layout = d3.layout.tree();
        var all_nodes = layout.nodes(root);
        var tval = all_nodes.map(function(d) {return +d[cUnit];});
        this.color = d3.scale.quantile()
            .domain([d3.min(tval), d3.max(tval)])
            .range(colors);
    };
    
    this.get_cmap = function(){
        var cc = this.color;
        var cs = this.colors;
        var m = this.colors.map(function (c){
            return ({color: c, value: d3.round(d3.mean(cc.invertExtent(c)), 2)});
        });
        return m;
    };

};

var LegendComponent = React.createClass({
    
    width: 200,
    height: 200,
    dispatcher: null,
    getInitialState : function (){
        return ({legend_created: false});
    }, 

    render : function(){
        return <svg ref="legend_svg" class="treelegend-container" id="treelegend_container" />
    },
    
    componentDidMount: function () {
        this.width = this.refs.svg.getDOMNode().offsetWidth;
        this.height = this.refs.svg.getDOMNode().offsetHeight;
    },
    
    componentWillUpdate : function(){
        if (!this.props.root) return false;
        var el = this.getDOMNode();
        if (this.state.legend_created){
            TreeLegend.update(el, this.props.appState);
        }else{
            this.dispatcher = TreeLegend.create(el, 
                {width: this.width, height:this.height}, 
                this.props.appState);
            this.state.legend_created = true;
        }
    }
});

///////////// Rendering the D3 tree
var Tree = React.createClass({
    
    dispatcher: null,
    width:200,
    height:200,
    getInitialState: function(){
        return ({
            tree_initialized:false,
            
        });
    },

    render: function() {
        return (
            <svg ref="svg" class="treeplot-container" id="treeplot_container">
            </svg>
        );
    },
    
    componentDidMount: function () {
        this.width = this.refs.svg.getDOMNode().offsetWidth;
        this.height = this.refs.svg.getDOMNode().offsetHeight;
        console.log(this.width)
    },

    componentDidUpdate : function(){

        console.log("Will update Tree view");
        if (!this.props.root) return false;
        
        var el = this.getDOMNode();
        if (this.state.tree_initialized){
            PhyloTree.update(el, this.props.appState, this.dispatcher);
        }else{
            var width = this.getDOMNode().offsetWidth;
            var height = this.getDOMNode().offsetHeight;

            var dispatcher = PhyloTree.create(el, {
            root:this.props.root,
            width:width, height:height}, 
            this.props.appState);
            this.dispatcher = dispatcher;
            dispatcher.on('point:tip_mouseover', this.select_tip);
            dispatcher.on('point:tip_mouseout', this.unselect_tip);
            this.setState({tree_initialized:true});
        }
        return false;
    },

    select_tip : function(d){
        console.log("TIP selected")
        console.log(d);
        this.props.setAppState({selected_tip : d.name});
    }, 

    unselect_tip : function(d){
        this.props.setAppState({selected_tip : null});
    }, 
});

var TreeTime = React.createClass({

    handleCheck : function(){
        console.log("Treetime CB changed");
        var tt = this.props.appState.treetime
        var xUnit = (!tt) ? "numdate" : "xvalue"
        this.props.setAppState({xUnit:xUnit, treetime:!tt});
    },

    render :function(){
        return (
            <div>
                <input type="checkbox" 
                onChange={this.handleCheck}
                checked={this.props.appState.treetime}></input>
                Show time axis
            </div>
            );
    }
});

var Mu = React.createClass({
    
    dispatcher: null, 
    
    
    render : function(){
        return <svg ref="mu_svg" class="mu-container" id="mu_container" />
    },

    componentDidMount: function () {
    },

    getInitialState : function(){
        return ({mu_initialized:false});
    },

    componentDidUpdate : function(){
        
        if (!this.props.mu) return false;
        var width = this.getDOMNode().offsetWidth;
        var height = this.getDOMNode().offsetHeight;
        var el = this.getDOMNode();
        if (!this.state.mu_initialized){
            this.dispatcher = MuPlot.create(el, 
                {
                    width:width, 
                    height:height, 
                    mu:this.props.mu
                }, 
                this.props.appState);
            this.dispatcher.on('point:point_mouseover', this.select_point);
            this.dispatcher.on('point:point_mouseout', this.unselect_point);
            this.state.mu_initialized = true;
        }else{
            MuPlot.update(el, this.props.mu, this.props.appState, this.dispatcher)
        }
    
    },

    select_point : function(d){
        this.props.setAppState({selected_tip:d.name})
    },

    unselect_point : function(d){
        this.props.setAppState({selected_tip:null})
    }

});

var RootLh = React.createClass({
    dispatcher: null, 
    render: function(){
        return <svg ref="lrooth_svg" class="lrooth-container" id="rootlh_container" />

    }, 
    componentDidUpdate : function(){
        
        if (!this.props.lh) return false;
        var width = this.getDOMNode().offsetWidth;
        var height = this.getDOMNode().offsetHeight;
        var el = this.getDOMNode();
        if (!this.props.appState.root_lh_initialized){
            this.dispatcher = RootLhPlot.create(el, 
                {
                    width:width, 
                    height:height, 
                    lh:this.props.lh
                }, 
                this.props.appState);
            this.dispatcher.on('point:point_mouseover', this.select_point);
            this.dispatcher.on('point:point_mouseout', this.unselect_point);
            this.props.setAppState({root_lh_initialized : true});
        }else{
            RootLhPlot.update(el, this.props.lh, this.props.appState, this.dispatcher)
        }
    
    },

    select_point : function(d){
        
    },

    unselect_point : function(d){
        
    }

});

///// Main APP 
var Res = React.createClass({
    
    root:{},
    mu:null,
    lh:null, 
    getInitialState : function (){
        return ({
            treetime:false, 
            cscale: new CScale('numdate'),
            xUnit:'xvalue',
            selected_tip:null,
            root_lh_initialized :false,
        });
    },

    componentWillMount : function(){
        console.log("Will get the tree file now. ");
        d3.json("/sessions/" + this.props.UID + "/out_tree.json", this.on_root);
        d3.json("/sessions/" + this.props.UID + "/out_tips.json", this.on_mu);
        d3.json("/sessions/" + this.props.UID + "/out_root_lh.json", this.on_root_lh);
    },
    on_root : function (err, root){
        if (err){console.warn("Can not get root node for the tree"); return;}
        this.root = root;
        this.state.cscale.update(this.state.cscale.cUnit, this.root);
        this.forceUpdate()
    },

    on_root_lh : function (err, lh){
        if (err){console.warn("Can not get root node for the tree"); return;}
        this.lh = lh;
        this.forceUpdate()
    },
    
    on_mu : function(err, mu){
        this.mu=mu;
        this.forceUpdate();
    },

    select_tip : function(d){
        console.log("Tip selected: " + d.strain);
    },
    
    unselect_tip : function(d){
        console.log("Tip unselected: " + d.strain);
    },

    setAppState :function (partialState, callback){
        console.log(partialState);
        return this.setState(partialState, callback);
    },

    render : function(){
        return (
            <div id="main_container">
            <h1>Results</h1>
            <div id="left_col">
            <TreeTime 
                appState={this.state}
                setAppState={this.setAppState}/>

            <LegendComponent 
                  root={this.root}
                  appState={this.state}
                  setAppState={this.setAppState}/>

            </div> 

            <div id="right_col">
            <h3>Phylogenetic tree</h3>
            <Tree UID={this.props.UID} 
                  root={this.root}
                  appState={this.state}
                  setAppState={this.setAppState}/>
            
            <h3>Mutation rate calculation:</h3>
            <Mu   mu={this.mu}
                appState={this.state}
                setAppState={this.setAppState}/>
            <h3>Likelihood distribution for the Tmrca:</h3>
            <RootLh lh={this.lh}
                appState={this.state}
                setAppState={this.setAppState}/>
            </div>
            </div>
        );
    }, 
    on_treetime_changed : function(){
        var checked = this.state.treetime;
        console.log(checked)
        this.setState({treetime : !checked})
        console.log("Treetime Cb changed to " + this.state.treetime)
        this.setState({xUnit: this.state.treetime ? "tvalue" : "xvalue"});
    }
});




export default Res;