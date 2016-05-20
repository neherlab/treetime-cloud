import React from  'react'
import ReactDOM from 'react-dom'
var request = require('superagent');
import { Panel, Button, Grid, Row, Col, FormControl, Checkbox, Table } from "react-bootstrap";

import PhyloTree from './phylo_tree.js'
import TreeLegend from './tree_legend.js'
import MuPlot from './mu_plot.js'
import RootLhPlot from './root_lh.js'

var Globals = require('./globals.js')
var colors = Globals.colors;

import Header from './header.js'
import Footer from './footer.js'

var DefaultScale = function(){
    this.get_color = function(x){
        return "#808080";
    }
    this.get_cmap = function(){
        return ([{color: "#808080", value:"ALL"}])
    }
};

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

///////////// Rendering the D3 tree
var TreeContainer = React.createClass({
    componentWillUpdate : function(){
        if (!this.props.root){
            return false;
        }
    },

    render: function(){
        return (
            <div className="results-section-wide" id="results-section_tree">
                <TreeLeftPane
                    root={this.props.root}
                    appState={this.props.appState}
                    setAppState={this.props.setAppState}
                    resetLayout={this.resetLayout}/>
                <TreeRightPane className="results-section-right_pane" id="results-section_tree-right_pane"
                    ref="TreeRightPane"
                    root={this.props.root}
                    appState={this.props.appState}
                    setAppState={this.props.setAppState}/>
            </div>
        );
    },
    resetLayout:function(){
        console.log("RessetLayout")
        var trp = this.refs.TreeRightPane
        trp.resetLayout()
    },

});

var TreeLeftPane = React.createClass({

    getInitialState : function(){
        return ({
            pos_disabled : true,
            pos_selected: 1,

        });
    },

    render: function(){
        return (
            <div className="results-section-left_pane" id="results-section_tree-left_pane">
                <h2>Phylogenetic tree</h2>
                <TreeTime
                    root={this.props.root}
                    appState={this.props.appState}
                    setAppState={this.props.setAppState}/>

                <Button bsStyle="primary" onClick={this.props.resetLayout}>Reset Layout</Button>

                <h3>Color by: </h3>
                <div id="results-section_tree-left_pane-select_colorscheme" >
                <FormControl componentClass="select" placeholder="numdate_given" className="select-treetime" id="results-section_tree-left_pane-select_colorval"
                        onChange={this.scaleChanged}>
                    <option value="numdate_given">numdate_given</option>
                    <option value="N">Nucleotide</option>
                    {
                        this.props.appState.terminal_colorby.map(function(d) {
                            if (d!="numdate_given")return <option key={d} value={d}>{d}</option>;
                        })
                    }
                </FormControl>
                <div id="results-section_tree-left_pane-div_pos">
                    <span id="results-section_tree-left_pane-pos" >
                        Pos:
                    </span>

                        <FormControl type="number" className="treetime-input-number" id="results-section_tree-left_pane-select_pos"
                            onChange={this.posChanged}
                            disabled={this.state.pos_disabled}>
                    </FormControl>
                </div>

                </div>


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

        default:

            console.log(value)

            var cValFunc = function(d){
                var md = d.terminal_metadata.filter(function(d){return d.name==value})
                if (md.length == 0) return null;
                return md[0].value;
            }
            var tips = []
            PhyloTree.gatherTips(this.props.root, tips)
            var all_values = tips.map(cValFunc)

            // get all types of the metadata entries
            all_values = all_values.filter(function(d, index){return all_values.indexOf(d) == index})
            var all_types = all_values.map(function(d){return typeof(d)})
            all_types = all_types.filter(function(d, index){return all_types.indexOf(d) == index;})

            var cScale;
            if (all_types.length == 1 && all_types[0] == "number"){ // all numbers - continuous scale
                cScale = new CScale();

            }else{ // strings or messed types - unique color for ever value
                cScale = new CategorialScale();

            }

            cScale.create(all_values);
            this.props.setAppState({
                    cvalue : cValFunc,
                    cscale: cScale
                })

            this.setState({pos_disabled : true});
            break;
        }
        this.forceUpdate();
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

    componentDidUpdate : function(){
        if (!this.props.root) return false;
        var el = this.getDOMNode();

        if (this.state.legend_created){
            TreeLegend.update(el, this.props.appState);
        }else{
            this.dispatcher = TreeLegend.create(el,
                {style:this.styleLegend},
                this.props.appState);
            this.state.legend_created = true;
        }
        return true;
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
                <Checkbox
                onChange={this.handleCheck}
                checked={this.props.appState.treetime}>
                Toggle time-tree
                </Checkbox>
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

    resetLayout : function(){
        if (!this.state.tree_initialized){return;}
        var el = this.getDOMNode();
        PhyloTree.resetLayout(el, this.dispatcher);
    },

    render: function() {
        return (
                <div
                    id="results-section_tree-right_pane"
                    className="results-section-right_pane"
                    ref="tree_svg"/>
        );
    },

    componentDidMount: function () {

    },

    componentDidUpdate : function(){

        //console.log("Will update Tree view");
        if (!this.props.root) {
            console.log("Results page cannot update Tree Container: No tree root defined")
            return false;
        }

        var el = this.getDOMNode();
        if (this.state.tree_initialized){
            PhyloTree.update(el, this.props.appState, this.dispatcher);
        }else{

            this.setState({tree_initialized:true});

            var dispatcher = PhyloTree.create(el, {
                root:this.props.root},
                this.props.appState);

            this.dispatcher = dispatcher;
            if (!dispatcher || typeof dispatcher == 'undefined') return;

            dispatcher.on('point:tip_mouseover', this.select_tip);
            dispatcher.on('point:tip_mouseout', this.unselect_tip);
        }
        return true;
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
            <div className="results-section results-section-narrow" id="results-section_mu">
            {/*<MuLeftPane
                appState={this.props.appState}
                setAppState={this.props.setAppState}
                mu={this.props.mu}
                root={this.props.root}/>*/}
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
            <div className="results-section-left_pane" id="results-section_mu-left_pane">
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
        return <div
            className="results-section results-section-narrow" id="results-section_mu"
            ref="mu_svg"/>
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
            <div  className="results-section results-section-narrow" id="results-section_tmrca">
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
            <div className="results-section-left_pane" id="results-section_tmrca-left_pane">
            <h2>LH distribution for tree root</h2>
            </div>
            );
    }
});

var TmrcaRightPane = React.createClass({

    dispatcher: null,
    render: function(){
        return <div
            className="results-section results-section-narrow" id="results-section_tmrca"
            ref="lrooth_svg"/>

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

var DownloadContainer = React.createClass({
    render: function(){
        return(
            <div className="results-section results-section-wide" id="results-section_dowload">
                <div className="spacer"></div>
                <DownloadLeftPane
                    appState={this.props.appState}
                    setAppState={this.props.setAppState}/>

                <DownloadRightPane
                    appState={this.props.appState}
                    setAppState={this.props.setAppState}/>
            </div>
        )
    }
});

var DownloadLeftPane = React.createClass({
    render : function(){
        return(
            <div className="results-section-left_pane" id="results-section_download-left_pane">
                <a className="btn btn-primary btn-file" id="results-section_download-btn_download" href={"/sessions/" + this.props.appState.UID + "/treetime_results.zip"} target="_blank">Download results (.zip)</a>
            </div>
        )
    }

});

var DownloadRightPane = React.createClass({
    render : function(){
        return(
            <div className="results-section-right_pane" id="results-section_download-right_pane">

                <Panel collapsible defaultCollapsed header="List of files in download archive" className="panel-treetime" id="results-panel_download_files">
                    <Table>
                        <thead>
                            <th>File name</th>
                            <th>Description</th>
                        </thead>
                        <tbody>
                            <tr>
                                <td>out_newick_tree.nwk</td>
                                <td>Phylogenetic tree in newick format. The branches lengths in the tree  are optimized with the TreeTime algorithm. Each internal node is assigned a name in format "NODE_XXXXXX" to link the nodes to the additional data provided (see below). Due to the limitations of the newick format, there is no additional information available in thetree file.</td>
                            </tr>
                            <tr>
                                <td>out_aln.fasta</td>
                                <td>Multiple sequence alignment including the inferred sequences of the internal nodes of the tree.</td>
                            </tr>
                            <tr>
                                <td>out_metadata.csv</td>
                                <td>Extended metadata in csv format. This file contains the metadata provided by the user. On top of this, there are the inferred dates for all nodes of the tree ("numdate"), including the leaves. The dates are in the format of "YYYY.F". There is also the values for the relaxed molecular clock ("mutation_rate/avg"), the deviation of the branch lengs from the optimal values ("branch_len/opt") and time since MRCA in numeric date format.</td>
                            </tr>
                            <tr>
                                <td>out_gtr.csv</td>
                                <td>Matrix for GTR model of evolution, used  in the TreeTime run. If "Infer from tree" was chosen, than the GTR matrix returned is one inferred from the data. For the detailed description of how the model was inferred from the tree, see the methods page. See our GTR project documentation (coming soon) for more details.</td>
                            </tr>
                            <tr>
                                <td>out_root_lh.csv</td>
                                <td>Likelihood distribution of the tree root position. The distribution is normalized to one.</td>

                            </tr>
                            <tr>
                                <td>out_molecular_clock.csv</td>
                                <td>Molecular clock estimate from the input tree</td>
                            </tr>
                            <tr>
                                <td>settings.json</td>
                                <td>Settings file used to run TreeTime</td>
                            </tr>
                            <tr>
                                <td>out_tree.json</td>
                                <td>Phylogenetic tree in json format. The tree includes entire information of the alignment and metadata. It also stores some useful technical information for tree representation</td>
                            </tr>
                        </tbody>
                    </Table>

                </Panel>
            </div>
        )
    }

});
///// Main APP
var Results = React.createClass({

    root:undefined,
    mu:null,
    lh:null,

    getInitialState : function (){
        return ({
            treetime:false,
            cvalue : function(d){
                return "";
            },
            cscale: new DefaultScale(),
            xUnit:'xvalue',
            selected_tip:null,
            root_lh_initialized :false,
            color_nuc_pos: 1,
            terminal_colorby: [],
            internal_metadata: [],
        });
    },

    on_root : function (err, root){

        //console.log("ROOT node came")
        console.log(root)
        if (err){
            console.warn("Can not get root node for the tree");
            console.warn(err);
            return;
        }
        this.root = root;
        this._update_lh_from_root();


        // load data to color terminal nodes and branches
        var tips = []
        PhyloTree.gatherTips(root, tips)
        var all_metas = tips.map(function(t){
            return t.terminal_metadata.map(function (m){return m.name})
        })
        var merged = Array.from(new Set([].concat.apply([], all_metas)));
        this.setState({terminal_colorby:merged})

        // create initial legend
        var cValFunc = function(d){
            var md = d.terminal_metadata.filter(function(d){return d.name=="numdate_given"})
            if (md.length == 0) return null;
            return md[0].value;
        }
        var tips = []
        PhyloTree.gatherTips(root, tips)
        var all_values = tips.map(cValFunc)
        var cScale = new CScale(); // legend with continuous scale
        cScale.create(all_values);
        this.setState({
            cvalue : cValFunc,
            cscale: cScale
        })

        this.forceUpdate()

    },

    on_root_lh : function (err, lh){
        if (err){
            console.warn("Can not get root node for the tree");
            console.warn(err);
            return;
        }
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
        window.addEventListener("resize", this.setAppState);
    },

    render : function(){
        return (

            <div>
                <Header />
                <div className="hugespacer"></div>
                <div className="page_container page_wide" >
                <TreeContainer
                    UID={this.props.UID}
                    root={this.root}
                    appState={this.state}
                    setAppState={this.setAppState}/>

                <MuRightPane
                    mu={this.mu}
                    root={this.root}
                    appState={this.state}
                    setAppState={this.setAppState}/>

                <TmrcaRightPane
                    lh={this.lh}
                    appState={this.state}
                    setAppState={this.setAppState}/>
                <div className="hugespacer"></div>
                <DownloadContainer
                    appState={this.state}
                    setAppState={this.setAppState}/>
                </div>
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