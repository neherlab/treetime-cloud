import React from  'react'
import ReactDOM from 'react-dom'
import Collapsible from 'react-collapsible'
import { Panel, PanelGroup, Button, Grid, Row, Col, FormControl, Checkbox, Table } from "react-bootstrap";
var request = require('superagent');

import Header from './header.js'
import Footer from './footer.js'





var DoBuildTree = React.createClass({
    getInitialState(){
        return ({
            checked:this.props.AppConfig.do_build_tree
        }
        );
    },

    handleChange(e){
        var build = this.state.checked;
        this.state.checked = !build;
        //console.log("Build tree Checkbox state changed to: " + this.state.checked);
        this.props.SetAppConfig({do_build_tree:this.state.checked})
    },
    
    render(){
        return (
            <div>
              <Checkbox className="cbox-treetime" id="cbox-treetime-build_tree" 
                checked={this.props.AppConfig.do_build_tree}
                onChange={this.handleChange}> 
            Build new tree
            </Checkbox>
            </div>
        );
    } 
});

var ShouldReuseBranchLength = React.createClass({

    getInitialState(){
        return (
        {
            checked: this.props.AppConfig.reuse_branch_len
        }
        );
    },    
    
    handleChange(e){
        var chk = this.props.AppConfig.reuse_branch_len;
        this.state.checked = !chk;
        this.props.SetAppConfig({"reuse_branch_len" : !this.props.AppConfig.reuse_branch_len});
    },
    

    render(){
        return (
            <div className="config-entry">
               <Checkbox className="cbox-treetime" id="cbox-treetime-reuse_branches" 
                checked={this.props.AppConfig.reuse_branch_len}
                onChange={this.handleChange}>
                Reuse tree branches
                </Checkbox>
            </div>
        );
    }
});

var DoReRoot = React.createClass({

    
    getInitialState(){
        return (
        {
            checked:this.props.AppConfig.reroot
        }
        );
    },    
    
    handleChange(e){
        var chk = this.props.AppConfig.reroot
        this.state.checked = !chk;
        //console.log("doReRoot Checkbox state changed to: " + this.state.checked);
        this.props.SetAppConfig({"reroot":this.state.checked});
    },
    

    render(){
        return (
            <div className="config-entry">
               <Checkbox className="cbox-treetime" id="welcome-panel_config-cbox_reroot" 
                checked={this.props.AppConfig.reroot}
                onChange={this.handleChange}>
                Optimize tree root position
                </Checkbox>
            </div>
        );
    }

});

var GTRmodel = React.createClass({
    
    onChange: function(value){
        console.log(value.target.value)
        var val = value.target.value
        if (val == "InferFromTree"){
            this.props.SetAppConfig({"gtr":"jukes_cantor", "infer_gtr":true})
        }else{
            this.props.SetAppConfig({"gtr":val, "infer_gtr":false})
        }
    }, 

    render: function(){

        return (
            <div className="config-entry">
            <span className="treetime-span-input" id="welcome-panel_config-span_GTR">GTR model:</span>
            <FormControl componentClass="select" placeholder="InferFromTree" className="select-treetime" id="welcome-panel_config-select_GTR"
            onChange={this.onChange}>
              <option value= "InferFromTree">Infer from tree</option>
              <option value= "jukes_cantor">Jukes-Cantor</option>
            </FormControl>
            </div>
        );
    }
});

var UseBranchPenalty = React.createClass({

    getInitialState(){
        return (
        {
            checked: this.props.AppConfig.use
        }
        );
    },    

    validate(text){
        return true;
    },
      
    handleCBChange(){
        var bool = this.state.settings.bool;
        var settings = this.state.settings;
        settings.bool = !bool;
        this.state.settings = settings;
        //console.log("Use Branch len penalty changed to: " + settings.bool);
        //console.log(this.state.settings);
        this.props.settings.change_handle(this.props.settings.name, this.state.settings);
    },

    handleTextChange(e){
        //console.log(e.target.value);
        var text = e.target.value;
        this.validate(text);
        var settings = this.state.settings;
        settings.value = text;
        settings.bool = true;
        this.state.settings = settings;
        this.props.settings.change_handle(this.props.settings.name, this.state.settings);
    },

    render(){

        return (
            <div className="config-entry">

                <Checkbox className="cbox-treetime" id="welcome-panel_config-cbox_penalty"
                    checked={this.state.settings.bool}
                    onChange={this.handleCBChange}> 
                Branch penalty = 
                </Checkbox>
                
                

                <FormControl type="number" className="treetime-input-number" id="welcome-panel_config-input_penalty"
                    disabled={!this.state.settings.bool} 
                    onChange={this.handleTextChange}/> 

              
            </div>
        );
    }
});


var UseSlope = React.createClass({

    getInitialState: function(){
        return (
        {
            use_mu: this.props.AppConfig.use_mu,
            mu: this.props.AppConfig.mu
        }
        );
    },    

    validate(text){
        return true;
    },
      
    handleCBChange(){
        var checked = this.props.AppConfig.use_mu
        this.props.SetAppConfig({use_mu: !checked})
        if (!this.props.AppConfig.use_mu){
            this.props.SetAppConfig({mu: 0.0})
        }
    },

    handleTextChange(e){
        //console.log(e.target.value);
        var text = e.target.value;
        if (!this.validate(text)){
            text=0.0
        }
        this.props.SetAppConfig({mu:text})
    },

    render(){

        return (
            <div className="config-entry">
              <Checkbox className="cbox-treetime" id="welcome_panel_config-cbox_slope"
                    checked={this.props.AppConfig.use_mu}
                    onChange={this.handleCBChange}> 
                
                Mutation rate (&mu;) = 
              
              </Checkbox>
              <div className="div-block">
              <FormControl type="number" className="treetime-input-number div-block"
                    disabled={!this.props.AppConfig.use_mu} 
                    onChange={this.handleTextChange}
                    value={this.props.AppConfig.mu}/> 
              
              <span className="treetime-span-input">(#/year)</span>
              </div>
            
            </div>
        );
    }
});

var DoResolvePoly = React.createClass({

    getInitialState(){
        return (
        {
            checked: this.props.AppConfig.resolve_poly
        }
        );
    },    
    
    handleChange(e){
        var chk = this.props.AppConfig.resolve_poly
        this.state.checked = !chk;
        //console.log("DoResolvePoly Checkbox state changed to: " + this.state.checked);
        this.props.SetAppConfig({"resolve_poly":!this.props.AppConfig.resolve_poly})
    },
    

    render(){
        return (
            <div className="config-entry">
               <Checkbox className="cbox-treetime"
                checked={this.props.AppConfig.resolve_poly}
                onChange={this.handleChange}>
                    Resolve polytomies
                </Checkbox>
            </div>
        );
    }
});

var DoCoalescent = React.createClass({
    
    getInitialState(){
        return (
        {
            Tc: this.props.AppConfig.Tc,
            coalescent:this.props.AppConfig.coalescent
        }
        );
    },    

    validate(text){
        return true;
    },
      
    handleCBChange(){
        var  chk = this.props.AppConfig.coalescent
        this.props.SetAppConfig({"coalescent": !chk})
        if (chk){
            this.props.SetAppConfig({"Tc": 0.0})
        }
    },

    handleTextChange(e){
        //console.log(e.target.value);
        var text = e.target.value;
        if (!this.validate(text)){
            text = 0.0
        }
        this.props.SetAppConfig({"Tc":text});
    },

    render(){

        return (
            <div className="config-entry" id="welcome_panel_config-div_coalescent">
                
                    <Checkbox className="cbox-treetime" id="welcome_panel_config-cbox_coalescent"
                           checked={this.props.AppConfig.coalescent}
                           onChange={this.handleCBChange}>
        
                        Model coalescent process. 
                    
                    </Checkbox>
                
                    <div id="div-block">
                        
                        <span className="treetime-span-input"> Tc = </span> 
                        
                        <FormControl  type="number" className="treetime-input-number"
                            disabled={!this.props.AppConfig.coalescent} 
                            onChange={this.handleTextChange}
                            value={this.props.AppConfig.Tc}/> 
                        
                        <span  className="treetime-span-input" id="welcome_panel_config-coalescence_spanHammDist"> (Hamming distance) </span> 
                    
                    </div>
            
            </div>
        );
    }

});

var DoRelaxedClock = React.createClass({

    getInitialState(){
        return (
        {
            relax_mu: this.props.AppConfig.relax_mu,
            slack: this.props.AppConfig.slack,
            coupling: this.props.AppConfig.coupling,
        }
        );
    },    

    validate(text){
        return true;
    },
      
    handleCBChange(){
        var chk = this.props.AppConfig.relax_mu
        if (chk){
            this.props.SetAppConfig({
                relax_mu:!chk,
                coupling:0.0,
                slack:0.0})
        }else{
            this.props.SetAppConfig({relax_mu:!chk})
        }
    },

    handleBChange(e){
        //console.log(e.target.value);
        var text = e.target.value;
        if (!this.validate(text)){
            text = 0.0
        }
        this.props.SetAppConfig({"slack":text})
    },

    handleAChange(e){
        //console.log(e.target.value);
        var text = e.target.value;
        if (!this.validate(text)){
            text = 0.0
        }
        this.props.SetAppConfig({"coupling":text})

    },
       
    render(){

        return (
            <div className="config-entry"  id="welcome_panel_config-div_relaxed" >
              
                <Checkbox className="cbox-treetime" id="welcome_panel_config-cbox_relaxed"
                checked={this.props.AppConfig.relax_mu}
                onChange={this.handleCBChange}>
                    Relax molecular clock
                </Checkbox>
              
                <div className="div-block">
                
                <span className="treetime-span-input" id="welcome_panel_config-beta_relaxed"> 
                    Slack(&alpha;) = 
                </span> 

                 <FormControl className="treetime-input-number" id="welcome_panel_config-bval_relaxed"
                        type= "text" 
                        disabled={!this.props.AppConfig.relax_mu}
                        onChange={this.handleBChange}
                        value={this.props.AppConfig.slack}/> 
                </div>

                <div className="div-block">
                
                    <span className="treetime-span-input"> 
                        Coupling(&beta;) = 
                    </span> 
                 
                    <FormControl className="treetime-input-number"
                        type= "number" 
                        disabled={!this.props.AppConfig.relax_mu}
                        onChange={this.handleAChange}
                        value={this.props.AppConfig.coupling}/>
                </div>

              
            </div>
        );
    }
});

var DoRootJoint = React.createClass ({
   
    getInitialState(){
        return (
        {
            checked: this.props.settings.settings[this.props.settings.name]
        }
        );
    },    
    
    handleChange(e){
        var chk = this.state.checked;
        this.state.checked = !chk;
        //console.log("doRootJoint Checkbox state changed to: " + this.state.checked);
        this.props.settings.change_handle(this.props.settings.name, this.state.checked);
    },
    

    render(){
        return (
            <div className="config-entry">
               <Checkbox className="cbox-treetime" id="welcome_panel_config_variance"
                checked={this.state.checked}
                onChange={this.handleChange}>
                Compute Root variance
                </Checkbox>
            </div>
        );
    }
});


//#TODO other properties (define them!)

var TreeTimeForm = React.createClass({

    setAppState : function(partialState){
        this.setState(partialState);
    },
    
    getInitialState : function() {
      return {
        UID: "JHG", 
        config: {
            do_build_tree:false,
            reuse_branch_len:false,
            reroot:false, 
            use_mu:false,
            coalescent:false,
            relax_mu:false,
            mu:0.0
        },
        tree_file:false,
        tree_filename:"Select tree file",
        aln_file:false,
        aln_filename:"Select alignment file",
        meta_file:false,
        meta_filename:"Select meta data file",
      };
    },

    handle_run: function(){
        //console.log("APP:: RUN button pressed");
        if ((!this.state.tree_file & !this.state.config.do_build_tree) || ! this.state.aln_file || !this.state.meta_file){
          var msg = "Cannot proceed with TreeTime: one or more file not loaded.\n\n"
          if ((!this.state.tree_file & !this.state.config.do_build_tree)){
            msg += "Phylogenetic tree file is missing.\n\n";
          }
          if (!this.state.aln_file){
            msg += "Sequence alignment file is missing.\n\n";
          }
          if (!this.state.meta_file){
            msg += "Meta data file is missing.\n\n";
          }
          alert(msg);
          return;
        }
        request.post("/" + this.state.UID)
          .set('Content-Type', 'application/json')
          .send({config: this.state.config})
          .end(this.on_run_status);
    },

    on_run_status : function(err, res){
    
        //console.log("RUN RESPONSE");
        //console.log(res)
        if (err){
            alert("Cannot start TreeTime calculations. Server error.")
            console.log(err)
            return;
        }
        window.location.replace("/" + this.state.UID + "/progress");

    },

    on_settings_changed : function(name, setting){
      //console.log("APP:: settings changed. " + name + " new value = " + setting);
      var settings = this.state.config
      settings[name] = setting;

      
      //this.state.settings = settings;

      if (name == "do_build_tree" ){
        if (setting){
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
      if (name=="GTR"){
            console.log(setting)
      }
      console.log(this.state.config)
    },


    componentDidMount: function(){
        //console.log("Welcome has been mounted");
        var parentNode = this.getDOMNode().parentNode;
        var UID = (parentNode.attributes.userid.value);
        //console.log("UID: " + UID)
        this.state.UID = UID;
        //console.log(this.state);
        console.log(cfg)
        this.SetAppConfig(cfg)
    },

    validate_form : function(){
        //console.log("Validating the form...");
        return null;
    },
    
    uploadTreeFile :function(evt){
        
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
        request.post('/upload/' + this.state.UID + '/file')
            .send(formData)
            .end(this.on_upload_tree);
    },

    on_upload_tree: function(err, res){
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
        request.post('/upload/' + this.state.UID + '/file')
            .send(formData)
            .end(this.on_upload_aln);
    },
    on_upload_aln: function(err, res){
        if (err){
            this.setAppState({aln_file:false});
            alert("Alignment file upload error. Please, try once more.")
            return;
        }
        this.setState({aln_filename:JSON.parse(res.text).AlnFile, aln_file:true})

    },

    uploadMetaFile :function(evt){

        if (evt.target.files.length == 0){
            // console.log("Resetting treefile")
            // this.setState({tree_filename:"No file chosen"})
            return;
        }
        //console.log("Uploading metadata file...");
        var formData = new FormData();
        formData.append('metafile', evt.target.files[0]);
        //for (var key in evt.target.files) {
        //    formData.append(key, files[key]);
        //}

        this.setAppState({meta_file:true});
        request.post('/upload/' + this.state.UID + '/file')
            .send(formData)
            .end(this.on_upload_meta);
    },
    on_upload_meta: function(err, res){
        if (err){
            this.setAppState({meta_file:false});
            alert("Meta data file upload error. Please, try once more.")
            return;
        }

        this.setState({meta_filename:JSON.parse(res.text).MetaFile, meta_file:true})

    },

    on_example_H3N2_NA_20 : function(){
        this.runExample("H3N2_NA_20");
    },

    on_example_H3N2_NA_500 : function(){
        this.runExample("H3N2_NA_500");
    },

    runExample : function(example){
        console.log("run example requested" + example)
        request.post("/" + this.state.UID + "/example")
          .set('Content-Type', 'application/json')
          .send({example:example})
          .end(this.on_example_upload_status);
    },

    on_example_upload_status : function(err, res){
    
        console.log("Example files uploaded")
        this.on_upload_aln(err, res)
        this.on_upload_meta(err, res)
        this.on_upload_tree(err, res)
    
    }, 

    SetAppConfig : function(cfg){
        console.log(cfg)
        var new_config = this.state.config
        for (var key in cfg) { 
            new_config[key] = cfg[key]; 
            if (key == "do_build_tree" ){
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
        this.setState({config:new_config})
    console.log(this.state.config)
    },

    render:function(){
        //console.log(this.state.settings)
        return (
            <div>
                <Header/>
                <div className="page_container">
    
                    <h2>Run TreeTime on server</h2>

                        <Panel collapsible defaultExpanded header="Upload data" className="panel-treetime" id="welcome_panel_files">
                            <Grid id="welcome_upload_grid"> 

                            <Row className="grid-treetime-row"> 
                            
                                <Col  xs={6} md={4} 
                                    id="welcome_col_upload_tree" className="grid-treetime-col-right" >
                                    
                                    <span className="btn btn-primary btn-file btn-file-treetime" id="btn-1">
                                        Newick <input  type="file" disabled={this.state.config.do_build_tree}  onChange={this.uploadTreeFile} />
                                    </span>
                                    
                                    {this.state.tree_filename}
                                    
                                    <DoBuildTree
                                        AppConfig={this.state.config}
                                        SetAppConfig={this.SetAppConfig}
                                        
                                        />

                                </Col>
                            </Row> 
                            
                            <Row className="grid-treetime-row"> 
                                                                
                                <Col  xs={6} md={4} className="grid-treetime-col-right">
                                    
                                    <span className="btn btn-primary btn-file btn-treetime btn-file-treetime">
                                        Fasta <input type="file" onChange={this.uploadAlnFile} />
                                    </span>
                                    {this.state.aln_filename}
                                
                                </Col>
                            </Row> 

                            <Row className="grid-treetime-row"> 
                                
                                <Col xs={6} md={4} className="grid-treetime-col-right"> 
                                    <span className="btn btn-primary btn-file btn-treetime btn-file-treetime">
                                        CSV <input type="file" onChange={this.uploadMetaFile} />
                                    </span>
                                    {this.state.meta_filename}
                                </Col>
                            </Row> 
                            </Grid> 
                        </Panel>
    
                        <Panel collapsible defaultCollapsed header="Alternatively, choose example dataset" className="panel-treetime" id="welcome_panel_examples">
                            <Row>
                            <Col xs={6} md={4} className="grid-treetime-col-right">
                                <button className="link_button" onClick={this.on_example_H3N2_NA_20}>Influenza H3N2 NA, 20 sequences</button>
                            </Col> 
                            <Col xs={6} md={4} className="grid-treetime-col-right">
                            </Col> 
                            </Row> 

                            <Row> 
                            
                            <Col xs={6} md={4} className="grid-treetime-col-right">
                                <button className="link_button" onClick={this.on_example_H3N2_NA_500}>Influenza H3N2 NA, 500 sequences</button>
                            </Col> 
                            
                            <Col xs={6} md={4} className="grid-treetime-col-right">
                            </Col> 

                            </Row> 
                        </Panel>

                    <Panel collapsible defaultCollapsed header="Advanced configuration" className="panel-treetime" id="welcome_panel_config">
                                
                        <ShouldReuseBranchLength 
                            AppConfig={this.state.config}
                            SetAppConfig={this.SetAppConfig}/>
                        
                        <DoReRoot 
                            AppConfig={this.state.config}
                            SetAppConfig={this.SetAppConfig}/>
        
                        <GTRmodel 
                            AppConfig={this.state.config}
                            SetAppConfig={this.SetAppConfig}/>
        
        
                        <UseSlope  
                            AppConfig={this.state.config}
                            SetAppConfig={this.SetAppConfig}/>
        
                        <DoResolvePoly 
                            AppConfig={this.state.config}
                            SetAppConfig={this.SetAppConfig}/>
        
                        <DoCoalescent  
                            AppConfig={this.state.config}
                            SetAppConfig={this.SetAppConfig}/>
                        
                        <DoRelaxedClock  
                            AppConfig={this.state.config}
                            SetAppConfig={this.SetAppConfig}/>
        
                        
                    </Panel>

                    <Button bsStyle="primary" className="btn-treetime" onClick={this.handle_run}>
                    Run TreeTime
                    </Button> 


                <div id="welcome_description">
                    <p>Welcome to the TreeTime server.
                    The description and HOWTO will appear here shortly. Please scroll down in order to 
                    run tree-time on the server.</p>
                </div>
                </div>
                <Footer/>
            </div>
        );
    }
});



/////////////// rendering

ReactDOM.render((
    <TreeTimeForm/>),
    document.getElementById('react'));


export default TreeTimeForm;