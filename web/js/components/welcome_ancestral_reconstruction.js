import React from  'react'
import ReactDOM from 'react-dom'
import Header from './header.js'
import GTR from './gtr.js'
var request = require('superagent');
import { Panel, PanelGroup, Button, Grid, Row, Col, FormControl, FormGroup, ControlLabel , Checkbox, Table, OverlayTrigger, Tooltip } from "react-bootstrap";

var PanelFiles = React.createClass({

    onBuildTreeSelected: function(evt){
        var chk = this.props.TreeAncConfig.build_tree != false && this.props.TreeAncConfig.build_tree != null;
        var btree = !chk;
        this.props.setTreeAncConfig({"build_tree":btree})
    },

    getInitialState(){
        return ({
                available_gtrs:[]
            }
        );
    },

    componentWillUpdate(nextProps, nextState) {
        console.log("Will Update: ");
        console.log(nextProps)

        var new_gtrs = nextProps.TreeAncConfig.available_gtrs;
        if (new_gtrs.length != this.state.available_gtrs.length){
            this.setState ({available_gtrs:nextProps.TreeAncConfig.available_gtrs})
        }
    },

    render: function(){
        return (
            <div>
            <Panel collapsible defaultExpanded header="Upload data" className="panel-treetime" id="welcome_panel_files">
            <Grid id="welcome_upload_grid">

                {/*Upload Tree file*/}
                <Row className="grid-treetime-row">
                    <Col xs={6} md={4} id="welcome_col_upload_tree" className="grid-treetime-col-right" >
                        <span className="btn btn-primary btn-file btn-file-treetime" id="btn-1">
                            Newick
                            <input type="file" disabled={this.props.TreeAncConfig.build_tree} onChange={this.props.uploadTreeFile}/>
                        </span>
                        {this.props.appState.tree_filename}
                        <Checkbox
                            checked={this.props.TreeAncConfig.build_tree}
                            onChange={this.onBuildTreeSelected}>
                            Build tree
                        </Checkbox>
                    </Col>
                </Row>

                {/*Upload Fasta file*/}
                <Row className="grid-treetime-row">
                    <Col  xs={6} md={4} className="grid-treetime-col-right">
                        <span className="btn btn-primary btn-file btn-treetime btn-file-treetime">
                            Fasta
                            <input type="file" onChange={this.props.uploadAlnFile} />
                        </span>
                        {this.props.appState.aln_filename}
                    </Col>
                </Row>

                {/*GTR model*/}
                <Row className="grid-treetime-row">
                    <Col xs={6} md={4} className="grid-treetime-col-right">
                    <GTR AppState={this.state} setTreeAncConfig={this.props.setTreeAncConfig} setGtrState={this.props.setGtrState}/>
                    </Col>
                </Row>

            </Grid>
            </Panel>
            </div>
        );
    }
});

var WelcomeAncPage = React.createClass({

    componentDidMount: function(){
        var parentNode = this.getDOMNode().parentNode;
        // set user id from the HTML template
        this.setAppState({'UID':parentNode.attributes.userid.value});
        this.setState({'treeAncConfig':cfg})
    },

    setAppState : function(partialState){
        this.setState(partialState);
    },

    getInitialState : function() {
      return {
        UID: null,
        // labels and status of the files uploads
        tree_file:false,
        tree_filename:"Select tree file",
        aln_file:false,
        aln_filename:"Select alignment file",
        // treeanc configuration
        treeAncConfig: {}
      };
    },

    setTreeAncConfig: function(cfg){
        var new_config = this.state.treeAncConfig
        for (var key in cfg) {

            new_config[key] = cfg[key];
            if (key == "build_tree" ){
                if (cfg[key]){
                    this.setState({
                    tree_filename: "Will be built from alignment",
                    tree_file: false,
                })
            }else{
                this.setState({
                    tree_filename: "Select tree file",
                    tree_file: false,
                })
            }
        }
        }
        this.setState({treeAncConfig:new_config})
        console.log("TreeAnc config changed. New config: ");
        console.log(this.state.treeAncConfig)
    },

    // Files uploads section
    uploadTreeFile :function(evt){

        this.setTreeAncConfig({"build_tree":false})

        if (evt.target.files.length == 0){
            console.log("Resetting treefile")
            this.setState({
                tree_filename:"Select tree file",
                tree_file:false,
            });
            return;
        }
        console.log("Uploading tree file...");
        var formData = new FormData();
        formData.append('treefile', evt.target.files[0]);
        //for (var key in evt.target.files) {
        //    formData.append(key, files[key]);
        //}


        this.setAppState({tree_file:true});
        request.post('/upload/treetime/' + this.state.UID + '/file')
            .send(formData)
            .end(this.onUploadTreeFile);
    },

    onUploadTreeFile: function(err, res){
        if (err){
            this.setAppState({tree_file:false, tree_filename: "Error uploading file"});
            alert("Tree file upload error. Please, try once more.")
            return;
        }

        this.setState({
            tree_filename:JSON.parse(res.text).TreeFile,
            tree_file:true
        })
    },


    uploadAlnFile :function(evt){

        if (evt.target.files.length == 0){
            // console.log("Resetting treefile")
            // this.setState({tree_filename:"No file chosen"})
            return;
        }
        //console.log("Uploading alignment file...");
        var formData = new FormData();
        formData.append('alnfile', evt.target.files[0]);
        //for (var key in evt.target.files) {
        //    formData.append(key, files[key]);
        //}
        this.setAppState({aln_file:true});
        request.post('/upload/treetime/' + this.state.UID + '/file')
            .send(formData)
            .end(this.onUploadAlnFile);
    },

    onUploadAlnFile: function(err, res){
        if (err){
            this.setAppState({aln_file:false});
            alert("Alignment file upload error. Please, try once more.")
            return;
        }
        this.setState({aln_filename:JSON.parse(res.text).AlnFile, aln_file:true})

    },


    onRunTreeAnc: function(){
        //console.log("APP:: RUN button pressed");
        if ((!this.state.tree_file & !this.state.treeAncConfig.build_tree) || ! this.state.aln_file){
          var msg = "Cannot proceed with TreeAnc: one or more file not loaded.\n\n"
          if ((!this.state.tree_file & !this.state.treeAncConfig.build_tree)){
            msg += "Phylogenetic tree file is missing.\n\n";
          }
          if (!this.state.aln_file){
            msg += "Sequence alignment file is missing.\n\n";
          }
          alert(msg);
          return;
        }
        request.post("/ancestral/" + this.state.UID + "/run")
          .set('Content-Type', 'application/json')
          .send({config: this.state.treeAncConfig})
          .end(this.onRunTreeAncResponse);
    },

    onRunTreeAncResponse : function(err, res){

        //console.log("RUN RESPONSE");
        //console.log(res)
        if (err){
            alert("Cannot start TreeAnc calculations. Server error.")
            console.log(err)
            return;
        }
        window.location.replace("/ancestral/" + this.state.UID + "/progress");
    },

    setGtrState: function(key, param_name, param_value){
        console.log("Welcome page: setting GTR state")
        var cfg = this.state.treeAncConfig
        var gtr = cfg.available_gtrs[key]
        if (!gtr.params){
            alert("Cannot set GTR parameter: This GTR has no parameters.")
            return;
        }
        for (var i = 0; i < gtr.params.length; ++i){
            var param = gtr.params[i];
            if (param.name == param_name){
                this.state.treeAncConfig.available_gtrs[key].params[i].value = param_value
                this.forceUpdate()
                break;
            }
        }
    },

    render:function(){
        return (
            <div>
                <Header/>
                <div className="page_container">
                <PanelFiles
                    TreeAncConfig={this.state.treeAncConfig}
                    setTreeAncConfig={this.setTreeAncConfig}
                    setGtrState={this.setGtrState}
                    appState={this.state}
                    uploadTreeFile={this.uploadTreeFile}
                    uploadAlnFile={this.uploadAlnFile}/>

                    <Button bsStyle="primary" className="btn-treetime" onClick={this.onRunTreeAnc}>
                        Run ancestral reconstruction
                    </Button>
                </div>

            </div>
        );
    }
});

ReactDOM.render((
    <WelcomeAncPage/>),
    document.getElementById('react'));

export default WelcomeAncPage;
