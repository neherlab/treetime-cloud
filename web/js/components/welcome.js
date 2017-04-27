import React from  'react'
import ReactDOM from 'react-dom'
import { Panel, PanelGroup, Button, Grid, Row, Col, FormControl, Checkbox, Table, OverlayTrigger, Tooltip } from "react-bootstrap";
var request = require('superagent');
import Header from './header.js'
import Footer from './footer.js'

var DoBuildTree = React.createClass({
    getInitialState(){
        return ({
            checked:this.props.AppConfig.build_tree
        }
        );
    },

    handleChange(e){
        var build = this.props.AppConfig.build_tree
        this.state.checked = !build;
        //console.log("Build tree Checkbox state changed to: " + this.state.checked);
        this.props.SetAppConfig({build_tree: !this.props.AppConfig.build_tree})
    },



    render(){
        const tooltip = (
            <Tooltip id="tooltip">
               Select to build a tree using <strong>FastTree</strong>.
               If more sophisticated methods are needed, provide your own tree.
            </Tooltip>
        );

        return (
            <div>

              <Checkbox className="cbox-treetime" id="cbox-treetime-build_tree"
                checked={this.props.AppConfig.build_tree}
                onChange={this.handleChange}>
                <OverlayTrigger placement="bottom" overlay={tooltip}>
                    <div> Build tree from alignment </div>
                </OverlayTrigger>
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
        const tooltip = (
            <Tooltip id="tooltip">
                Set checked if you are sure that the branch lenghts of the input tree
                are consitent with the GTR model (e.g. represent Hamming distance, or similar).
            </Tooltip>
        );

        return (
            <div className="config-entry">
               <Checkbox className="cbox-treetime" id="cbox-treetime-reuse_branches"
                checked={this.props.AppConfig.reuse_branch_len}
                onChange={this.handleChange}>
                    <OverlayTrigger placement="top" overlay={tooltip}>
                        <div>Reuse tree branches</div>
                    </OverlayTrigger>
                </Checkbox>
            </div>
        );
    }
});

var DoReRoot = React.createClass({


    getInitialState(){
        return (
        {
            checked:this.props.AppConfig.reroot == 'best'
        }
        );
    },

    handleChange(e){
        var chk = this.props.AppConfig.reroot
        this.state.checked = !chk;
        //console.log("doReRoot Checkbox state changed to: " + this.state.checked);
        this.props.SetAppConfig({reroot:this.state.checked ? 'best' : null});
    },


    render(){
        const tooltip = (
            <Tooltip id="tooltip">
                Re-root tree to optimal root-to-tip regression.
            </Tooltip>
        );
        return (
            <div className="config-entry">
               <Checkbox className="cbox-treetime" id="welcome-panel_config-cbox_reroot"
                checked={this.props.AppConfig.reroot}
                onChange={this.handleChange}>
                    <OverlayTrigger placement="top" overlay={tooltip}>
                    <div>Optimize tree root position</div>
                    </OverlayTrigger>
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
            this.props.SetAppConfig({"gtr":"infer"})
        }else{
            this.props.SetAppConfig({"gtr":"Jukes_Cantor"})
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



                <FormControl type="text" className="treetime-input-number" id="welcome-panel_config-input_penalty"
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

        const tooltip = (
            <Tooltip id="tooltip">
                Use fixed mutation rate instead of inferred clock rate.
            </Tooltip>
        );

        return (
            <div className="config-entry">
              <Checkbox className="cbox-treetime" id="welcome_panel_config-cbox_slope"
                    checked={this.props.AppConfig.use_mu}
                    onChange={this.handleCBChange}>

                    <OverlayTrigger placement="top" overlay={tooltip}>
                        <div>Fix the mutation rate. &mu; =</div>
                    </OverlayTrigger>

              </Checkbox>
              <div className="div-block">
              <FormControl type="text" className="treetime-input-number div-block"
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
                    Resolve polytomies using temporal constraints
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
            coalescent:this.props.AppConfig.model_coalescent
        }
        );
    },

    validate(text){
        return true;
    },

    handleCBChange(){
        var  chk = this.props.AppConfig.model_coalescent
        this.props.SetAppConfig({model_coalescent: !chk})
        if (chk){
            this.props.SetAppConfig({Tc: 0.0})
        }
    },

    handleTextChange(e){
        //console.log(e.target.value);
        var text = e.target.value;
        if (!this.validate(text)){
            text = 0.0
        }
        this.props.SetAppConfig({Tc:text});
    },

    render(){

        return (
            <div className="config-entry" id="welcome_panel_config-div_coalescent">

                    <Checkbox className="cbox-treetime" id="welcome_panel_config-cbox_coalescent"
                           checked={this.props.AppConfig.model_coalescent}
                           onChange={this.handleCBChange}>

                        Use coalescent prior. T<sub>c</sub> =

                    </Checkbox>

                    <div id="div-block">

                        <span className="treetime-span-input" id="welcome_panel_config_tc" >  </span>

                        <FormControl  type="text" className="treetime-input-number"
                            disabled={!this.props.AppConfig.model_coalescent}
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
            do_relaxed_clock: this.props.AppConfig.do_relaxed_clock,
            slack: this.props.AppConfig.relax_clock.slack,
            coupling: this.props.AppConfig.relax_clock.coupling,
        }
        );
    },

    validate(text){
        return true;
    },

    handleCBChange(){
        var chk = this.props.AppConfig.do_relaxed_clock
        if (chk){
            this.props.SetAppConfig({
                do_relaxed_clock:!chk,
                coupling:0.0,
                slack:0.0})
        }else{
            this.props.SetAppConfig({do_relaxed_clock:!chk})
        }
    },

    handleBChange(e){
        //console.log(e.target.value);
        var text = e.target.value;
        if (!this.validate(text)){
            text = 0.0
        }
        this.setState({slack:text})
        this.props.SetAppConfig({relax_clock:
            {
                slack:this.state.slack,
                coupling:this.state.coupling
            }
        })
    },

    handleAChange(e){
        //console.log(e.target.value);
        var text = e.target.value;
        if (!this.validate(text)){
            text = 0.0
        }
        this.setState({coupling:text})
        this.props.SetAppConfig({relax_clock:
            {
                slack:this.state.slack,
                coupling:this.state.coupling
            }
        })

    },

    render(){

        const tooltip = (
            <Tooltip id="tooltip">
                Perform an autocorrelated relaxation of the mutation rate along the tree.
            </Tooltip>
        );


        return (
            <div className="config-entry"  id="welcome_panel_config-div_relaxed" >

                <Checkbox className="cbox-treetime" id="welcome_panel_config-cbox_relaxed"
                checked={this.props.AppConfig.do_relaxed_clock}
                onChange={this.handleCBChange}>
                    <OverlayTrigger placement="top" overlay={tooltip}>
                    <div>Relax molecular clock</div>
                    </OverlayTrigger>
                </Checkbox>

                <div className="div-block">

                <span className="treetime-span-input" id="welcome_panel_config-beta_relaxed">
                    Slack: &alpha; =
                </span>

                 <FormControl className="treetime-input-text" id="welcome_panel_config-bval_relaxed"
                        type= "text"
                        disabled={!this.props.AppConfig.do_relaxed_clock}
                        onChange={this.handleBChange}
                        value={this.props.AppConfig.slack}/>
                </div>

                <div className="div-block">

                    <span className="treetime-span-input">
                        Coupling: &beta; =
                    </span>

                    <FormControl className="treetime-input-number"
                        type= "number"
                        disabled={!this.props.AppConfig.do_relaxed_clock}
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
            build_tree:false,
            reuse_branch_len:false,
            reroot:'best',
            use_mu:false,
            model_coalescent:false,
            do_relaxed_clock:false,
            relax_clock:{slack:1, coupling:1},
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
        if ((!this.state.tree_file & !this.state.config.build_tree) || ! this.state.aln_file || !this.state.meta_file){
          var msg = "Cannot proceed with TreeTime: one or more file not loaded.\n\n"
          if ((!this.state.tree_file & !this.state.config.build_tree)){
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

        this.SetAppConfig({"build_tree":false})

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

    on_example_zika_65 : function(){
        this.runExample("zika_65");
    },

    runExample : function(example){
        console.log("run example requested:  " + example)
        this.SetAppConfig({"build_tree":false})
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
        this.setState({config:new_config})
    console.log(this.state.config)
    },

    render:function(){
        //console.log(this.state.settings)
        const csv_tooltip = (
            <Tooltip className="csv_tooltip">

                Upload comma-separated file with sampling dates and additional metadata.
                The format of the file should be as follows:

                    <div className='spacer'/>

                    <table>
                    <tr>
                        <th>
                        name,date,location,...
                        </th>
                    </tr>
                    <tr>
                        <td>
                        A/Oregon/15/2009,2009.4,Oregon,... <br/>
                        A/New_York/182/2000,2000.1,New_York,...
                        </td>
                    </tr>
                    </table>

                    <div className='spacer'/>


                <ul className="scriptFont">
                <li>A single header row with column names is required. "name" and "date" are required fields.</li>

                <li>The first column should contain the <strong>name</strong> of the sequences.</li>

                <li>The first column with "date" in the name is used as tip dates.
                Admissible formats are: (i) numeric date as "2015.7", (ii) date string, e.g. "YYYY-MM-DD".
                </li>

                <li>Other columns may contain arbitrary data of any format and have arbitrary names.
                The server will try to parse these columns and show it in the results section.</li>
                </ul>





            </Tooltip>
        );


        return (

            <div>
                <Header/>


                <div className="page_container">
                <div id="welcome_description">

                <h4>Features</h4>
                    <ul>
                        <li>Ancestral sequence reconstruction</li>
                        <li>Molecular clock tree inference</li>
                        <li>Inference of GTR models</li>
                        <li>Rerooting to obtain best root-to-tip regression</li>
                        <li>Auto-correlated relaxed molecular clock (with normal prior)</li>
                    </ul>

                TreeTime source code is available on <a href='https://github.com/neherlab/treetime'>Github</a>.


                </div>

                <div className="bigspacer"/>

                        <Panel collapsible defaultExpanded header="Upload data" className="panel-treetime" id="welcome_panel_files">
                            <Grid id="welcome_upload_grid">

                            <Row className="grid-treetime-row">

                                <Col  xs={6} md={4}
                                    id="welcome_col_upload_tree" className="grid-treetime-col-right" >

                                    <span className="btn btn-primary btn-file btn-file-treetime" id="btn-1">
                                        Newick <input  type="file"  onChange={this.uploadTreeFile} />
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
                                    <OverlayTrigger placement="bottom" overlay={csv_tooltip}>
                                    <span className="btn btn-primary btn-file btn-treetime btn-file-treetime">
                                            CSV <input type="file" onChange={this.uploadMetaFile} />
                                    </span>
                                    </OverlayTrigger>
                                    {this.state.meta_filename}

                                </Col>
                            </Row>
                            </Grid>
                        </Panel>

                        <Panel collapsible defaultCollapsed header="Example datasets" className="panel-treetime" id="welcome_panel_examples">
                            <Table striped condensed hove>
                            <thead>
                              <tr>
                                <th>Specie</th>
                                <th>Region</th>
                                <th>Seq Len</th>
                                <th>#Seq</th>
                                <th>Dates range</th>
                                <th>Load</th>
                              </tr>
                            </thead>
                            <tbody>
                            <tr>
                              <th>Influenza H3N2</th>
                              <th>NA</th>
                              <th>1409</th>
                              <th>20</th>
                              <th>2000-2013</th>
                              <th><Button bsStyle="primary" className="btn-treetime" onClick={this.on_example_H3N2_NA_20}>Load</Button></th>
                            </tr>
                            <tr className="info-treetime">
                              <th>Influenza H3N2</th>
                              <th>NA</th>
                              <th>1409</th>
                              <th>500</th>
                              <th>1968-2010</th>
                              <th><Button bsStyle="primary" className="btn-treetime" onClick={this.on_example_H3N2_NA_500}>Load</Button></th>
                            </tr>
                            <tr>
                              <th>Zika</th>
                              <th>Full genome</th>
                              <th>10617</th>
                              <th>65</th>
                              <th>2013-2016</th>
                              <th><Button bsStyle="primary" className="btn-treetime" onClick={this.on_example_zika_65}>Load</Button></th>
                            </tr>
                            </tbody>
                            </Table>

                        </Panel>




                    <Panel collapsible defaultCollapsed header="Advanced configuration" className="panel-treetime" id="welcome_panel_config">



                        <DoReRoot
                            AppConfig={this.state.config}
                            SetAppConfig={this.SetAppConfig}/>

                        <DoResolvePoly
                            AppConfig={this.state.config}
                            SetAppConfig={this.SetAppConfig}/>

                        <GTRmodel
                            AppConfig={this.state.config}
                            SetAppConfig={this.SetAppConfig}/>


                        <UseSlope
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
