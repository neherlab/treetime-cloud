import React from  'react'
var request = require('superagent');

var DoBuildTree = React.createClass({
    getInitialState(){
        return ({
            checked: false,//this.props.settings.settings[this.props.settings.name]
        }
        );
    },

    handleChange(e){
        var build = this.state.checked;
        this.state.checked = !build;
        console.log("Build tree Checkbox state changed to: " + this.state.checked);
        //this.props.settings.change_handle(this.props.settings.name, this.state.checked);
    },
    
    render(){
        return (
            <div>
              <input type= "checkbox"
                checked={this.state.checked}
                onChange={this.handleChange}/> 
            Build tree with FastTree
            </div>
        );
    } 
});

var ShouldReuseBranchLength = React.createClass({

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
        console.log("Reuse Branches Checkbox state changed to: " + this.state.checked);
        this.props.settings.change_handle(this.props.settings.name, this.state.checked);
    },
    

    render(){
        return (
            <div>
               <input type= "checkbox"
                checked={this.state.checked}
                onChange={this.handleChange}>
                </input>
                Reuse tree branches
            </div>
        );
    }
});

var DoReRoot = React.createClass({

    
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
        console.log("doReRoot Checkbox state changed to: " + this.state.checked);
        this.props.settings.change_handle(this.props.settings.name, this.state.checked);
    },
    

    render(){
        return (
            <div>
               <input type= "checkbox"
                checked={this.state.checked}
                onChange={this.handleChange}>
                </input>
                Optimize tree root position
            </div>
        );
    }

});

var GTRmodel = React.createClass({
    render(){

        return (
            <div class="select">
            GTR model:    <select value= "A" >
              <option value= "A">Jukes-Cantor</option>
            </select>
            </div>
        );
    }
});

var UseBranchPenalty = React.createClass({

    getInitialState(){
        return (
        {
            settings: this.props.settings.settings[this.props.settings.name]
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
        console.log("Use Branch len penalty changed to: " + settings.bool);
        console.log(this.state.settings);
        this.props.settings.change_handle(this.props.settings.name, this.state.settings);
    },

    handleTextChange(e){
        console.log(e.target.value);
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
            <div>
              <div style={{floating:'left'}}> 
                <input type="checkbox" name="penalty_cb"
                    checked={this.state.settings.bool}
                    onChange={this.handleCBChange}> </input>
                Branch penalty = 
                <input style={{'margin-left':'10px'}} type="text" name="penalty_value"
                    disabled={!this.state.settings.bool} 
                    onChange={this.handleTextChange}/> 
              </div>
              
              <div>
                
              </div>
            </div>
        );
    }
});


var UseSlope = React.createClass({

    getInitialState(){
        return (
        {
            settings: this.props.settings.settings[this.props.settings.name]
        }
        );
    },    

    validate(text){
        return true;
    },
      
    handleCBChange(){
        var bool = this.state.settings.bool;
        var settings = this.state.settings; // copy 
        settings.bool = !bool;
        this.state.settings = settings;
        console.log("Use Slpe changed to: " + settings.bool);
        console.log(this.state.settings);
        this.props.settings.change_handle(this.props.settings.name, this.state.settings);
    },

    handleTextChange(e){
        console.log(e.target.value);
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
            <div>
              <input type="checkbox" name="use_slope"
                checked={this.state.settings.bool}
                onChange={this.handleCBChange}
              > </input>
              Mutation rate (&mu;) = 
              <input style={{'margin-left':'10px', 'margin-right':'10px'}}  type="text" name="slope_value"
                    disabled={!this.state.settings.bool} 
                    onChange={this.handleTextChange}/> 
              (#/year)
            </div>
        );
    }
});

var DoResolvePoly = React.createClass({

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
        console.log("DoResolvePoly Checkbox state changed to: " + this.state.checked);
        this.props.settings.change_handle(this.props.settings.name, this.state.checked);
    },
    

    render(){
        return (
            <div>
               <input type= "checkbox"
                checked={this.state.checked}
                onChange={this.handleChange}>
                </input>
                Resolve polytomies
            </div>
        );
    }
});

var DoCoalescent = React.createClass({
    
    getInitialState(){
        return (
        {
            settings: this.props.settings.settings[this.props.settings.name]
        }
        );
    },    

    validate(text){
        return true;
    },
      
    handleCBChange(){
        var bool = this.state.settings.bool;
        var settings = this.state.settings; // copy 
        settings.bool = !bool;
        this.state.settings = settings;
        console.log("Do coalescent changed to: " + settings.bool);
        console.log(this.state.settings);
        this.props.settings.change_handle(this.props.settings.name, this.state.settings);
    },

    handleTextChange(e){
        console.log(e.target.value);
        var text = e.target.value;
        this.validate(text);
        var settings = this.state.settings;
        settings.Tc = text;
        settings.bool = true;
        this.state.settings = settings;
        this.props.settings.change_handle(this.props.settings.name, this.state.settings);
    },

    render(){

        return (
            <div id="welcome_coalescent">
                <div>
                    <input type="checkbox" name="do_coalescent"
                           checked={this.state.settings.bool}
                           onChange={this.handleCBChange}></input>
        
                    Model coalescent process. 
                </div>
                <div style={{'margin-left':'10px', 'margin-right':'10px'}}>
                    Tc = <input style={{'margin-left':'10px', 'margin-right':'10px'}} 
                          type="text" name="tc"
                          disabled={!this.state.settings.bool} 
                          onChange={this.handleTextChange}/> 
                    (Hamming distance)
                </div>
            
            </div>
        );
    }

});

var DoRelaxedClock = React.createClass({

    getInitialState(){
        return (
        {
            settings: this.props.settings.settings[this.props.settings.name]
        }
        );
    },    

    validate(text){
        return true;
    },
      
    handleCBChange(){
        var bool = this.state.settings.bool;
        var settings = this.state.settings; // copy 
        settings.bool = !bool;
        this.state.settings = settings;
        console.log("Do coalescent changed to: " + settings.bool);
        console.log(this.state.settings);
        this.props.settings.change_handle(this.props.settings.name, this.state.settings);
    },

    handleAChange(e){
        console.log(e.target.value);
        var text = e.target.value;
        this.validate(text);
        var settings = this.state.settings;
        settings.alpha = text;
        settings.bool = true;
        this.state.settings = settings;
        this.props.settings.change_handle(this.props.settings.name, this.state.settings);
    },

    handleBChange(e){
        console.log(e.target.value);
        var text = e.target.value;
        this.validate(text);
        var settings = this.state.settings;
        settings.beta = text;
        settings.bool = true;
        this.state.settings = settings;
        this.props.settings.change_handle(this.props.settings.name, this.state.settings);
    },
       
    render(){

        return (
            <div id="welcome_relaxed" >
              <div style={{'margin-right':'10px'}} >
                <input type= "checkbox"
                checked={this.state.settings.bool}
                onChange={this.handleCBChange}>
                </input>
                Estimate autocorrelated molecular clock
              </div>
              <div style={{'margin-left':'10px', 'margin-right':'10px'}}  >
                &alpha; = 
                 
                 <input style={{'margin-left':'10px', 'margin-right':'10px'}} 
                        type= "text" 
                        disabled={!this.state.settings.bool}
                        onChange={this.handleAChange}/>
              
                &beta; = 
                 <input style={{'margin-left':'10px', 'marginRight':'10px'}} 
                        type= "text" 
                        disabled={!this.state.settings.bool}
                        onChange={this.handleBChange}/> 
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
        console.log("doRootJoint Checkbox state changed to: " + this.state.checked);
        this.props.settings.change_handle(this.props.settings.name, this.state.checked);
    },
    

    render(){
        return (
            <div>
               <input type= "checkbox"
                checked={this.state.checked}
                onChange={this.handleChange}>
                </input>
                Compute Root variance
            </div>
        );
    }
});

var DoCalcGTR = React.createClass ({
   
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
        console.log("doCalcGTR Checkbox state changed to: " + this.state.checked);
        this.props.settings.change_handle(this.props.settings.name, this.state.checked);
    },
    

    render(){
        return (
            <div>
               <input type= "checkbox"
                checked={this.state.checked}
                onChange={this.handleChange}>
                </input>
                Calc GTR from tree
            </div>
        );
    }
});

//#TODO other properties (define them!)

var TreeTimeForm = React.createClass({
    
    componentDidMount(){
        console.log(this.props);
    },

    handle_run(){
        
        console.log("Welcome:: RunButton pressed");
        this.validate_form();
        this.props.handle_run();
    
    },

    validate_form(){
        console.log("Validating the form...");
        return null;
    },
    
    on_settings_changed(name, settings){
        console.log("Welcome:: Setings changed: " + name + ".  new value: " + settings);
        this.props.handle_settings_change(name, settings);
    },

    uploadTreeFile (evt){
    
        console.log("Uploading tree file...");
        var formData = new FormData();
        formData.append('treefile', evt.target.files[0]);
        //for (var key in evt.target.files) {
        //    formData.append(key, files[key]);
        //}

        request.post('/' + this.props.UID + '/tree_file')
            .send(formData)
            .end(function(err, res){
                console.log("Got upload response!");
            });
    },

    uploadAlnFile (evt){
    
        console.log("Uploading alignment file...");
        var formData = new FormData();
        formData.append('alnfile', evt.target.files[0]);
        //for (var key in evt.target.files) {
        //    formData.append(key, files[key]);
        //}

        request.post('/' + this.props.UID + '/aln_file')
            .send(formData)
            .end(function(err, res){
                console.log("Got upload response!");
            });
    },

    uploadMetaFile (evt){

        console.log("Uploading metadata file...");
        var formData = new FormData();
        formData.append('metafile', evt.target.files[0]);
        //for (var key in evt.target.files) {
        //    formData.append(key, files[key]);
        //}

        request.post('/' + this.props.UID + '/meta_file')
            .send(formData)
            .end(function(err, res){
                console.log("Got upload response!");
            });
    },

    render(){
        var act = "/" + this.props.UID + "/run";
        console.log(act)
        console.log(this.props.settings)
        return (
            <div id="welcome_container">
                <h2>Welcome</h2>
                <div id='welcome_welcome'>
                <p>Welcome to the TreeTime server.
                The description and HOWTO will appear here shortly. Please scroll down in order to 
                run tree-time on the server.</p>
                </div>
                <h2>Run TreeTime on server:</h2>
            
                <h3>1. Upload data</h3>

                <div id="welcome_files">
                    
                    <div id="welcome_treeupload">
                    <h4 class="welcome-reeupload-header"> Upload tree file: </h4>
                    <input id="welcome_input_tree"  
                        type="file" 
                        name="treefile"
                        //disabled={this.props.settings.doBuildTree}
                        onChange={this.uploadTreeFile}></input>

                    <h4 id='welcome_treeupload_or'> Or: </h4>
                    <DoBuildTree settings={{
                        name: "doBuildTree",
                        settings:this.props.settings,
                        change_handle: this.on_settings_changed
                    }}/>
                    </div> 
                    
                    <div class="welcome_treeupload_header" id="welcome_alnupload">
                    <h4 > Upload alignment file: </h4>
                    <input  id="welcome_input_aln"  
                        type="file" 
                        name="alnfile"
                        onChange={this.uploadAlnFile}></input>
                    </div>
                    
                    <div id="welcome_metaupload">
                    <h4 class="welcome-reeupload-header"> Upload metadata: </h4>
                    <input  id="welcome_input_meta"  
                        type="file" 
                        name="metafile"
                        onChange={this.uploadMetaFile}></input>
                    </div>
                </div>

                <h3>2. Configure parameters</h3>
            
                <div id="welcome_params">

                <ShouldReuseBranchLength settings={{
                    name: "shouldReuseBranchLen",
                    settings:this.props.settings,
                    change_handle: this.on_settings_changed
                
                }}/>
                
                <DoReRoot settings={{
                    name: "doReroot",
                    settings:this.props.settings,
                    change_handle: this.on_settings_changed
                
                }}/>

                <GTRmodel settings={{
                    name: "GTRmodel",
                    settings:this.props.settings,
                    change_handle: this.on_settings_changed
                
                }}/>

                <UseBranchPenalty settings={{
                    name: "shouldUseBranchLenPenalty",
                    settings:this.props.settings,
                    change_handle: this.on_settings_changed
                
                }}/>

                <UseSlope settings={{
                    name: "shouldUseSlope",
                    settings:this.props.settings,
                    change_handle: this.on_settings_changed
                }}/>

                <DoResolvePoly settings={{
                    name: "doResolvePoly",
                    settings:this.props.settings,
                    change_handle: this.on_settings_changed
                }}/>

                <DoCoalescent settings={{
                    name: "doCoalescent",
                    settings:this.props.settings,
                    change_handle: this.on_settings_changed
                }}/>
                
                <DoRelaxedClock settings={{
                    name: "doRelaxedClock",
                    settings:this.props.settings,
                    change_handle: this.on_settings_changed
                }}/>
                
                <DoRootJoint settings={{
                    name: "doRootJoint",
                    settings:this.props.settings,
                    change_handle: this.on_settings_changed
                }}/>

                <DoCalcGTR settings={{
                    name: "doCalcGTR",
                    settings:this.props.settings,
                    change_handle: this.on_settings_changed
                }}/>

                </div>

                <div >
                    <input type='button' id="welcome_run" onClick={this.handle_run} />
                </div>

            </div>
        );
    }
});

export default TreeTimeForm;