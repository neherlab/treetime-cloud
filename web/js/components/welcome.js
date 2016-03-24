import React from  'react'

var DoBuildTree = React.createClass({
    getInitialState(){
        return ({
            checked: this.props.settings.settings[this.props.settings.name]
        }
        );
    },

    handleChange(e){
        var build = this.state.checked;
        this.state.checked = !build;
        console.log("Build tree Checkbox state changed to: " + this.state.checked);
        this.props.settings.change_handle(this.props.settings.name, this.state.checked);
    },
    
    render(){
        return (
            <div>
              <input type= "checkbox"
                checked={this.state.checked}
                onChange={this.handleChange}/> 
            Build tree
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
            <div>
            GTR model<select value= "A">
              <option value= "A">Jukes-Cantor</option>
            </select>
            </div>
        );
    }
});

var UseBranchPenalty = React.createClass({

    getInitialState(){
        return {use_penalty: true}
    },
  
    
    handleChange(){

        console.log("Penalty Checkbox state changed");
        this.setState({use_penalty : !this.state.use_penalty});
        return null;
    },


    render(){

        return (
            <div>
              <input type= "checkbox"
                checked={this.state.use_penalty}
                onChange={this.handleChange}
              > </input>
              Branch penalty
              <br/>
              <input type= "text" disabled={!this.state.use_penalty}/> Penalty value
            </div>
        );
    }
});


var UseSlope = React.createClass({

    getInitialState (){
        return {use_slope: true}
    },
        
    handleChange (){
        console.log("Slope Checkbox state changed");
        this.setState({use_slope : !this.state.use_slope});
        return null;
    },

    render(){
        return (
            <div>
              <input type= "checkbox"
                checked={this.state.should_build}
                onChange={this.handleChange}></input> 
              Slope date-mutation rate conversion
              <br/>
              <input type= "text" disabled={this.state.use_slope}/> Slope value (#muts/year)
            </div>
        );
    }
});

var DoResolvePoly = React.createClass({


    
    getInitialState (){
        return {use_penalty: true}
    },
  
    handleChange(){
        console.log("Penalty Checkbox state changed");
        this.setState({use_penalty : !this.state.use_penalty});
        return null;
    },
 

    render(){
        return (
            <div>
               <input type= "checkbox"
                checked={this.state.should_build}
                onChange={this.handleChange}>
                </input>
                Resolve polytomies
            </div>
        );
    }
});

var DoCoalescent = React.createClass({
    
    getInitialState(){
        return {do_coalescent: false}
    },
      
    handleChange(){
        console.log("Penalty Checkbox state changed");
        this.setState({do_coalescent : !this.state.do_coalescent});
        return null;
    },

    render(){
        return (
            <div>
              <input type= "checkbox"
                checked={this.state.do_coalescent}
                onChange={this.handleChange}>
                </input>
                Use coalescent model
                <br/>
                <input type= "text" disabled={!this.state.do_coalescent}/> Tc

            </div>
        );
    }

});

var DoRelaxedClock = React.createClass({

    
    getInitialState (){
        return {do_clock: false}
    },
   
    
    handleChange (){

        console.log("MolClock Checkbox state changed");
        this.setState({do_clock : !this.state.do_clock});
        return null;
    },
       
    render(){

        return (
            <div>
              <input type= "checkbox"
                checked={this.state.do_clock}
                onChange={this.handleChange}>
                </input>
                Estimate autocorrelated mutation rate 
                <br/>
                 <input type= "text" disabled={!this.state.do_clock}/> Alpha
                 <input type= "text" disabled={!this.state.do_clock}/> Beta
            </div>
        );
    }
});

var DoRootJoint = React.createClass ({

    
    getInitialState(){
        return {use_penalty: true}
    },
    
    
    handleChange(){

        console.log("RootJoint Checkbox state changed");
        this.setState({use_penalty : !this.state.use_penalty});
        return null;
    },
    

    render(){

        return (
            <div>
              <input type= "checkbox"
                checked={this.state.should_build}
                onChange={this.handleChange}></input> Compute joint root-slope distribution
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

    render(){
        var act = "/" + this.props.UID;
        console.log(act)
        console.log(this.props.settings)
        return (
            <div>
            <form  method='post' action={act} encType='multipart/form-data' >
            <div id="files">
                <h2>Upload files</h2>
                <input type="file" 
                       name="treefile"
                       disabled={this.props.settings.doBuildTree}></input>
                <input type="file" name="aln_file"></input>
                <input type="file" name="metafile"></input>
            </div>

            <div id="params">
                <DoBuildTree settings={{
                    name: "doBuildTree",
                    settings:this.props.settings,
                    change_handle: this.on_settings_changed
                }}/>

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

            </div>
            <input type='button' value="RUN treetime" onClick={this.handle_run} />
            </form>
            </div>
        );
    }
});

export default TreeTimeForm;