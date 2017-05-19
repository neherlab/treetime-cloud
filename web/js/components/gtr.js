import React from 'react'
import {FormControl, ControlLabel, Grid, Row, Col, OverlayTrigger, Tooltip} from "react-bootstrap";

var GTR = React.createClass({

    setAppState : function(partialState){
        this.setState(partialState, function(){
            console.log("GTR: state has been updated. Available gtrs: ", this.state.available_gtrs)
            console.log("GTR: state has been updated. Selected gtr: ", this.state.selected_gtr)
        });
    },

    getInitialState(){
        return ({
                available_gtrs:{},
                selected_gtr:null
            }
        );
    },

    componentDidUpdate(nextProps, nextState) {
        var new_gtrs = nextProps.AppState.available_gtrs;
        console.log("GTR panel will Update: ");
        var new_len = Object.keys(new_gtrs).length
        var old_len = Object.keys(this.state.available_gtrs).length
        if (new_len != old_len){
            this.setAppState ({available_gtrs:new_gtrs})
        }
    },

    onGtrSelected : function(e){
        var gtr = e.target.value
        console.log(this.state.available_gtrs[gtr])
        this.setAppState({selected_gtr: gtr})
        this.props.setTreeAncConfig({"gtr": gtr})
    },

    gtr_dropdown_params: function(key){
        if (! (key in this.state.available_gtrs) || key=='infer') {
            return null;
        }else{
            var value = this.state.available_gtrs[key].value;
            return <option key={key} value={key}>{value}</option>;
        }
    },

    //We cannot mutate elementDOM on-the-fly hence have to display all the params
    //of all GTR models and control their visibility
    gtrParamsVisibility: function(gtr){
        if (this.state.selected_gtr == gtr){
            return {"margin-top":"10px", "float":"none", "display":"block"}
        }else{
            return {"display":"none"}
        }
    },

    onParamChanged: function(key, param, evt){
        this.props.setGtrState(key, param, evt.target.value)
    },

    renderGtrParams: function(key){
        var gtr = this.state.available_gtrs[key];
        if (!gtr || !gtr.params){
            return (<div></div>);
        }else{
        var params  = []

        for (var i = 0; i < gtr.params.length; i++){
            params.push(this.renderGtrParam(key, gtr.params[i]));
        }
        return (
            <div style={this.gtrParamsVisibility(key)}>
            {
                params
            }
            </div>
        );
        }

    },

    renderGtrParam: function(key, param){
        var change_callback = function(evt){
            this.onParamChanged(key, param.name, evt)
        }.bind(this)

        var tip = param.tip
        var tooltip = (
        <Tooltip id="tooltip">
            {tip}
        </Tooltip>
        );
        var overlay={tooltip}

        var fc = (
            <FormControl
                type="number"
                step="1e-3"
                maxlength="5"
                min="0"
                max="1"
                disabled={false}
                style={{"display":"inline-block"}}
                onChange={change_callback}
                value={param.value}>

            </FormControl>
        );

        // in case there is no tip- do not show the tooltip
        var fc_component = tip ? (
            <OverlayTrigger placement="top"  overlay={tooltip}>
                <FormControl
                type="number"
                step="1e-3"
                maxlength="5"
                min="0"
                max="1"
                disabled={false}
                style={{"display":"inline-block"}}
                onChange={change_callback}
                value={param.value}>

            </FormControl>
            </OverlayTrigger>): fc;

        return (
            <span style={{"display":"inline-block", "float":"left", "margin-right":"10px", "margin-top":"10px"}} overlay={overlay}>

            <span style={{"display":"inline-block"}}>{param.name} =

            </span>
            {fc_component}
            </span>
        );
    },

    render: function(){
        return (
        <div>
            <ControlLabel>GTR model</ControlLabel>
            <FormControl componentClass="select"
                    placeholder="Infer from tree"
                    className="select-treetime"
                    id="welcome-panel_config-select_GTR"
                    onChange={this.onGtrSelected}>
                    <option key='infer' value='infer'>Infer from tree</option>
                {
                    Object.keys(this.state.available_gtrs).map(this.gtr_dropdown_params)
                }

                {
                // <option value= "infer">Infer from tree</option>
                // {
                //     this.state.available_gtrs.map(function(d){
                //         return <option key={d.key} value={d.key}>{d.value}</option>;
                //     })
                // }
                }
            </FormControl>

            {
                Object.keys(this.props.AppState.available_gtrs).map(this.renderGtrParams)
            }

        </div>
        )
    }
});

export default GTR;
