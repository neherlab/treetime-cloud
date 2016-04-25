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

var CScale = function(){

    this.colors = ["#4D92BF", "#5AA5A8", "#6BB18D", "#80B974", "#98BD5E", "#B1BD4E",
          "#C8B944", "#DAAC3D", "#E59738", "#E67732", "#E14F2A", "#DB2522"];
    
    this.create = function (all_values){
        this.color = d3.scale.quantile()
            .domain([d3.min(all_values), d3.max(all_values)])
            .range(colors);
    };

    this.get_color = function(x){
        if (!this.color){
            return "#808080";
        }else{
            return this.color(x);
        }
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

var NucScale = function(){
    
    this.colors = ["#4D92BF", "#98BD5E", "#E59738", "#DB2522"];
    this.nucs = ["A", "C", "G", "T"]
    
    this.get_color = function(val){

        //var val = color_val_getter(d);
        //console.log("NucScale generating color for value: " + val)
        var idx = this.nucs.indexOf(val)
        //console.log(idx)
        if (idx == -1){return "#808080";}
        else{ 
            //console.log(this.colors[idx])
            return this.colors[idx];
        }
    };

    this.get_cmap = function(){

        var m = [{color:this.colors[0], value: this.nucs[0]},
                 {color:this.colors[1], value: this.nucs[1]},
                 {color:this.colors[2], value: this.nucs[2]},
                 {color:this.colors[3], value: this.nucs[3]}]
        return m;
    
    };

};

var CategorialScale = function(){
    
    this.domain = [];
    this.cScale = null;

    this.create = function(domain){
        if (domain.length < 11){
            this.cScale = new d3.scale.category10()
            .domain(domain)
        }else{
            this.cScale = new d3.scale.category20b()
            .domain(domain)
        }
        this.domain = domain;
    }

    this.get_color = function(x){
        if (!this.cScale || this.domain.length == 0){
            return "#808080";
        }else{
            console.log(this.cScale)
            return this.cScale(x);
        }
    };
    
    this.get_cmap = function(){
        var get_color = this.cScale;
        return this.domain.map(function (c){
            return ({color: get_color(c), value: c});
        });
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
  'position':'relative',
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
    
    getInitialState : function(){
        return ({
            pos_disabled : false,
            pos_selected: 1,
        });
    }, 

    render: function(){
        return (
            <div style={styleLeftPane}>
                <h2>Phylogenetic tree</h2>
                <TreeTime 
                    root={this.props.root}
                    appState={this.props.appState}
                    setAppState={this.props.setAppState}/>
                <h3>Color by: </h3>
                
                <select onChange={this.scaleChanged}>
                    <option value="N">Nucleotide</option>
                    <option value="T">Sampling date</option>
                    <option value="C">Country</option>
                    <option value="R">Region</option>
                </select>

                <input type="number" style={{"margin-left":"30px"}}
                    onChange={this.posChanged}
                    disabled={this.state.pos_disabled}> Pos</input>

                <h3>Color codes:</h3>
                <LegendComponent 
                    root={this.props.root}
                    appState={this.props.appState}
                    setAppState={this.props.setAppState}/>
            </div>
            );
    }, 

    scaleChanged : function(value){
        var value = (value.target.value);
        switch(value){
        case ("T"):
            this.setState({pos_disabled : true});
            var cValFunc = function(d){return d.tvalue};
            var cScale = new CScale();
            
            var root = this.props.root
            if (typeof(root) == 'undefined') {
            
                // do nothing
            
            }else{
            
                var layout = d3.layout.tree();
                var all_nodes = layout.nodes(root);
                var tval = all_nodes.map(function(d) {return +cValFunc(d);});
                cScale.create(tval);
            
            }

            this.props.setAppState({
                    cvalue : cValFunc,
                    cscale: cScale
                })

            break;
        case ("N"):
            this.setState({pos_disabled : false});
            var cValFunc = function(d){return d.strseq[1]};
            //var cScale = new NucScale();
            var cScale = new CategorialScale();
            var values = ["A", "C", "G", "T"]
            values = values.filter(function(d, index){return values.indexOf(d) == index})
            cScale.create(values);
            this.props.setAppState({
                    cvalue : cValFunc,
                    cscale: cScale
                })

            break;
        case ("C"):
            this.setState({pos_disabled : true});
            var cValFunc = function(d){return d.country};
            var cScale = new CategorialScale();

            // create scale from all defined values 
            var root = this.props.root
            if (typeof(root) == 'undefined') {
                // do nothing
            }else{
                var layout = d3.layout.tree();
                var all_nodes = layout.nodes(root);
                var values = all_nodes.map(function(d) {return +cValFunc(d);});
                values = values.filter(function(d, index){return values.indexOf(d) == index})
                cScale.create(values);
            }

            this.props.setAppState({
                    cvalue : cValFunc,
                    cscale: cScale
                })
            break;
        case ("R"):
             this.setState({pos_disabled : true});
            var cValFunc = function(d){return d.region};
            var cScale = new CategorialScale();

            // create scale from all defined values 
            var root = this.props.root
            if (typeof(root) == 'undefined') {
                // do nothing
            }else{
                var layout = d3.layout.tree();
                var all_nodes = layout.nodes(root);
                var values = all_nodes.map(function(d) {return +cValFunc(d);});
                values = values.filter(function(d, index){return values.indexOf(d) == index})
                cScale.create(values);
            }

            this.props.setAppState({
                    cvalue : cValFunc,
                    cscale: cScale
                })
            break;
        default:
            this.setState({pos_disabled : true});
            break;
        }
    },

    posChanged: function(e){
        var pos_value = e.target.value;
        this.props.setAppState({
            cvalue : function(d){
                if (typeof(d) == 'undefined'){return;}
                return d.strseq[parseInt(pos_value, 10)];}
        });
    },

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
        //console.log("createing legend...");
        //console.log(this.props.root);
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
        //console.log("Treetime CB changed");
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

        //console.log("Will update Tree view");
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
            cvalue : function(d){
                return d.strseq[57];
            },
            cscale: new NucScale(),
            xUnit:'xvalue',
            selected_tip:null,
            root_lh_initialized :false,
            color_nuc_pos: 1,
        });
    },

    on_root : function (err, root){
        //console.log("ROOT node came")
          if (err){console.warn("Can not get root node for the tree"); return;}
        this.root = root;
        this._update_lh_from_root();
        this.forceUpdate()
    },

    on_root_lh : function (err, lh){
        if (err){console.warn("Can not get root node for the tree"); return;}
        this.lh = lh;
        //console.log("LH has been read from the server file...")
        this.forceUpdate()
    },

    _update_lh_from_root : function() {
        
    },
    
    on_mu : function(err, mu){
        this.mu=mu;
        this.forceUpdate();
    },

    select_tip : function(d){
        //console.log("Tip selected: " + d.strain);
    },
    
    unselect_tip : function(d){
        //console.log("Tip unselected: " + d.strain);
    },

    setAppState :function (partialState, callback){
        this.setState(partialState, callback);
        this.forceUpdate();
    },

    componentDidMount: function(){
        var parentNode = this.getDOMNode().parentNode;
        var UID = (parentNode.attributes.userid.value);
        this.state.UID = UID;
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
        this.setState({treetime : !checked})
        this.setState({xUnit: this.state.treetime ? "tvalue" : "xvalue"});
    }
});


ReactDOM.render((
    <Results />),
    document.getElementById('react'));



export default Results;