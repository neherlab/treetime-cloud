import React from  'react'
import ReactDOM from 'react-dom'
var request = require('superagent');
import PhyloTree from './phylo_tree.js'
import TreeLegend from './tree_legend.js'
import MuPlot from './mu_plot.js'
import RootLhPlot from './root_lh.js'

var Globals = require('./globals.js')
var colors = Globals.colors;

import Header from './header.js'
import Footer from './footer.js'

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

var styleMain = {
    float:'left',
    //background: "#eebbff",
    padding: "5px",
    margin: "5px",
    width:"100%",
    //height: "600px",
    //'border-style':'solid',
    //'border-width': '1px'
};

var styleLeftPane = {

  float:'left',
  width:'29.5%',

};

var styleRightPane = {
  height: "700px",
  //width:'400px',
  float:'left',
  width:'70%',
  //height:styleMain.height,
  'border-style':'solid',
  'border-width': '1px',
  //background: "#bbff66",
};

///////////// Rendering the D3 tree
var TreeContainer = React.createClass({
    render: function(){
        return (
            <div style={styleMain}>
                <TreeLeftPane 
                    root={this.props.root}
                    appState={this.props.appState}
                    setAppState={this.props.setAppState}/>
                <TreeRightPane 
                    root={this.props.root}
                    appState={this.props.appState}
                    setAppState={this.props.setAppState}/>
            </div>
        );
    }
});

var TreeLeftPane = React.createClass({
    
    render: function(){
        return (
            <div style={styleLeftPane}>
                <h2>Phylogenetic tree</h2>
                <TreeTime 
                    root={this.props.root}
                    appState={this.props.appState}
                    setAppState={this.props.setAppState}/>
                <h3>Color codes:</h3>
                <LegendComponent 
                    root={this.props.root}
                    appState={this.props.appState}
                    setAppState={this.props.setAppState}/>
            </div>
            );

    }

});

var LegendComponent = React.createClass({
    
    styleLegend :{
        width: "100%",
        height:"400px",
        position: "relative",
        //'border-style':'solid',
        //'border-width':'1px',
        rect_size: 30,
    },

    dispatcher: null,
    getInitialState : function (){
        return ({legend_created: false});
    }, 

    render : function(){
        return (
            <div  style={this.styleLegend} ref="legend_svg" class="treelegend-container" id="treelegend_container" />
            );
    },
    
    componentDidMount: function () {
    },
    
    componentWillUpdate : function(){
        console.log("createing legend...");
        console.log(this.props.root);
        if (!this.props.root) return false;
        var el = React.findDOMNode(this);
        
        if (this.state.legend_created){
            TreeLegend.update(el, this.props.appState);
        }else{
            this.dispatcher = TreeLegend.create(el, 
                {style:this.styleLegend}, 
                this.props.appState);
            this.state.legend_created = true;
        }
        return false; // prevent react from re-draw
    }
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
                Toggle time-tree
            </div>
            );
    }
});

var TreeRightPane = React.createClass({
    
    dispatcher: null,
    getInitialState: function(){
        return ({
            tree_initialized:false,
            
        });
    },

    render: function() {
        return (
            <div>
                <div ref="tree_svg" class="treeplot-container" id="treeplot_container" style={styleRightPane}/>
            </div>
        );
    },

    componentDidMount: function () {
        
    },

    componentDidUpdate : function(){

        console.log("Will update Tree view");
        if (!this.props.root) return false;
        
        var el = this.refs.tree_svg.getDOMNode();
        if (this.state.tree_initialized){
            PhyloTree.update(el, this.props.appState, this.dispatcher);
        }else{
            var width =  el.offsetWidth;
            var height = el.offsetHeight;

            var dispatcher = PhyloTree.create(el, {
                root:this.props.root}, 
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


var MuContainer = React.createClass({
    render: function(){
        return (
            <div style={styleMain}>
            <MuLeftPane 
                appState={this.props.appState}
                setAppState={this.props.setAppState}
                mu={this.props.mu}
                root={this.props.root}/>
            <MuRightPane 
                appState={this.props.appState}
                setAppState={this.props.setAppState}
                root={this.props.root}/>
            </div>
            );
    }
});

var MuLeftPane = React.createClass({
    
    mu : function(){
        return this.props.appState.mol_clock ? 
            this.props.appState.mol_clock.slope.toExponential(3) 
            : "---";
    }, 

    r2 : function(){
        return this.props.appState.mol_clock ? 
            Math.round(this.props.appState.mol_clock.r2 * 1000) / 1000
            : "---";
    }, 

    render: function(){
        return (
            <div style={styleLeftPane}> 
            <h2>Molecular clock</h2>
            <h4>Average mutation rate:<br/> &mu; = {this.mu()} year<sup>-1</sup></h4>
            <h4>Correlation coefficient:<br/> R<sup>2</sup> = {this.r2()}</h4>
            </div>
            );
    }
});

var MuRightPane = React.createClass({
    
    dispatcher: null, 
    
    
    render : function(){
        return <div style={styleRightPane} ref="mu_svg" class="mu-container" id="mu_container" />
    },

    componentDidMount: function () {
    },

    getInitialState : function(){
        return ({mu_initialized:false});
    },

    componentDidUpdate : function(){
        
        if (!this.props.root) return false;
        var el = this.getDOMNode();
        if (!this.state.mu_initialized){
            this.dispatcher = MuPlot.create(el, 
                {
                    root:this.props.root
                }, 
                this.props.appState);
            this.dispatcher.on('point:point_mouseover', this.select_point);
            this.dispatcher.on('point:point_mouseout', this.unselect_point);
            this.dispatcher.on('mol_clock:regression_changed', this.on_regression_changed)
            this.state.mu_initialized = true;
        }else{
            MuPlot.update(el, this.props.root, this.props.appState, this.dispatcher)
        }
    
    },

    on_regression_changed : function(regression){
        this.props.setAppState({mol_clock: regression});
        
    },

    select_point : function(d){
        this.props.setAppState({selected_tip:d.name})
    },

    unselect_point : function(d){
        this.props.setAppState({selected_tip:null})
    }

});


var TmrcaContainer = React.createClass({
    render: function(){
        return (
            <div  style={styleMain}>
            <TmrcaLeftPane
                lh={this.props.lh}
                appState={this.props.appState}
                setAppState={this.props.setAppState}/>
            <TmrcaRightPane
                lh={this.props.lh}
                appState={this.props.appState}
                setAppState={this.props.setAppState}/>
            </div>
            );
    }
});

var TmrcaLeftPane = React.createClass({
    render: function(){
        return (
            <div style={styleLeftPane}>
            <h2>Tmrca LH distribution</h2>
            </div>
            );
    }
});

var TmrcaRightPane = React.createClass({
    
    dispatcher: null, 
    render: function(){
        return <div style={styleRightPane} ref="lrooth_svg" class="lrooth-container" id="rootlh_container" />

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
var Results = React.createClass({
    
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

    on_root : function (err, root){
        console.log("ROOT node came")
          if (err){console.warn("Can not get root node for the tree"); return;}
        this.root = root;
        this._update_lh_from_root();
        this.state.cscale.update(this.state.cscale.cUnit, this.root);
        this.forceUpdate()
    },

    on_root_lh : function (err, lh){
        if (err){console.warn("Can not get root node for the tree"); return;}
        this.lh = lh;
        console.log("LH has been read from the server file...")
        this.forceUpdate()
    },

    _update_lh_from_root : function() {
        
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
    componentDidMount: function(){
        console.log("Welcome has been mounted");
        var parentNode = this.getDOMNode().parentNode;
        var UID = (parentNode.attributes.userid.value);
        console.log("UID: " + UID)
        this.state.UID = UID;
        console.log(this.state);
        console.log("Will get the tree file now. ");
        d3.json("/sessions/" + this.state.UID + "/out_tree.json", this.on_root);
        d3.json("/sessions/" + this.state.UID + "/out_tips.json", this.on_mu);
        d3.json("/sessions/" + this.state.UID + "/out_root_lh.json", this.on_root_lh);
        
    },

    render : function(){
        return (
            
            <div>
                <Header />
            
                <h1>Results</h1>
                <TreeContainer
                    UID={this.props.UID} 
                    root={this.root}
                    appState={this.state}
                    setAppState={this.setAppState}/>
               
                <MuContainer
                    mu={this.mu}
                    root={this.root}
                    appState={this.state}
                    setAppState={this.setAppState}/>
                
                <TmrcaContainer 
                    lh={this.lh}
                    appState={this.state}
                    setAppState={this.setAppState}/>
                <Footer/>
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


ReactDOM.render((
    <Results />),
    document.getElementById('react'));



export default Results;