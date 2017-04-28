import React from  'react'
import ReactDOM from 'react-dom'
import Header from './header.js'
var request = require('superagent');
import { Panel, PanelGroup, Button, Grid, Row, Col, FormControl, FormGroup, ControlLabel , Checkbox, Table, OverlayTrigger, Tooltip } from "react-bootstrap";

var PanelText = React.createClass({
    render: function(){
        return (
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


            </div>
        );
    }
});

var PanelFiles = React.createClass({
    render: function(){
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
                <Panel collapsible defaultExpanded header="Upload data" className="panel-treetime" id="welcome_panel_files">
                            <Grid id="welcome_upload_grid">

                            <Row className="grid-treetime-row">

                                <Col  xs={6} md={4}
                                    id="welcome_col_upload_tree" className="grid-treetime-col-right" >

                                    <span className="btn btn-primary btn-file btn-file-treetime" id="btn-1">
                                        Newick <input  type="file"  onChange={this.uploadTreeFile} />
                                    </span>

                                    TreeFilename{/*this.state.tree_filename*/}

                                    {// <DoBuildTree
                                    //     AppConfig={this.state.config}
                                    //     SetAppConfig={this.SetAppConfig}
                                    //     blah={'blah'}
                                    //     />
                                    }

                                </Col>
                            </Row>

                            <Row className="grid-treetime-row">

                                <Col  xs={6} md={4} className="grid-treetime-col-right">

                                    <span className="btn btn-primary btn-file btn-treetime btn-file-treetime">
                                        Fasta <input type="file" onChange={this.uploadAlnFile} />
                                    </span>
                                    AlnFilename{/*this.state.aln_filename*/}

                                </Col>
                            </Row>

                            <Row className="grid-treetime-row">

                                <Col xs={6} md={4} className="grid-treetime-col-right">
                                    <OverlayTrigger placement="bottom" overlay={csv_tooltip}>
                                    <span className="btn btn-primary btn-file btn-treetime btn-file-treetime">
                                            CSV <input type="file" onChange={this.uploadMetaFile} />
                                    </span>
                                    </OverlayTrigger>
                                    MetaFileName{/*this.state.meta_filename*/}

                                </Col>
                            </Row>
                            </Grid>
                        </Panel>
            </div>
        );
    }
});


var PanelExamples = React.createClass({
    render: function(){
        return (
            <div>
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

            </div>
        );
    }

});

var PanelConfig = React.createClass({
    render: function(){
        const reroot_tooltip = (
            <Tooltip id="tooltip">
                Re-root tree to optimal root-to-tip regression.
            </Tooltip>
        );

        return (

            <div>
                    <Panel collapsible defaultCollapsed header="Advanced configuration" className="panel-treetime" id="welcome_panel_config">
                        {/*Reroot to best root*/}
                        <Checkbox>
                            {/*checked={this.props.AppConfig.reroot}
                            onChange={this.handleChange}>*/}
                            <OverlayTrigger placement="top" overlay={reroot_tooltip}>
                            <div>Optimize tree root position</div>
                            </OverlayTrigger>
                        </Checkbox>

                        {/*Polytomies resolution*/}
                        <Checkbox>
                            {/*checked={this.props.AppConfig.resolve_poly}
                            onChange={this.handleChange}>*/}
                            Resolve polytomies using temporal constraints
                        </Checkbox>

                        {/*GTR model*/}
                        <FormGroup>
                            <ControlLabel>GTR model</ControlLabel>
                            <FormControl componentClass="select" placeholder="InferFromTree" className="select-treetime" id="welcome-panel_config-select_GTR"
                                onChange={this.onChange}>
                                <option value= "InferFromTree">Infer from tree</option>
                                <option value= "jukes_cantor">Jukes-Cantor</option>
                            </FormControl>
                        </FormGroup>
                        <FormGroup>
                            <Checkbox>Fix substitution rate</Checkbox>
                            <FormControl type="number"
                                disabled={false}
                                /*
                                onChange={this.handleTextChange}*/
                                value={0.001}/>
                            <span>(#/year)</span>
                        </FormGroup>

                        <FormGroup>
                            <Checkbox>Use coalescent prior</Checkbox>
                            <span style={{"visibility":"hidden", "display":"inline-block"}}>Tc = </span>
                            <FormControl style={{"visibility":"hidden", "display":"inline-block"}} type="number"
                                /*disabled=true
                                onChange={this.handleTextChange}*/
                                value={0.001}/>

                        </FormGroup>

                        <FormGroup>
                            <Checkbox>Relax molecular clock</Checkbox>
                            <div style={{"verticalAlign":"center", "lineHeight":"50px", "width":"100%","display":"inline-block", "float":"left"}}>
                            <span  style={{"display":"inline-block", "float":"left"}}>Slack: &alpha; =</span>
                            <FormControl type="number"  style={{"height":"50px", "display":"inline-block", "float":"left"}}
                                visible={false}/*disabled=true
                                onChange={this.handleTextChange}*/
                                value={0.02}/>
                            <span  style={{"display":"inline-block", "float":"left"}}>Coupling: &beta; =</span>
                            <FormControl type="number" style={{"visibility":"hidden", "display":"inline-block", "float":"left"}}
                                /*disabled=true
                                onChange={this.handleTextChange}*/
                                value={0.0001}/>
                            </div>
                        </FormGroup>


                        {/*}
                        <UseSlope
                            AppConfig={this.state.config}
                            SetAppConfig={this.SetAppConfig}/>

                        <DoCoalescent
                            AppConfig={this.state.config}
                            SetAppConfig={this.SetAppConfig}/>

                        <DoRelaxedClock
                            AppConfig={this.state.config}
                            SetAppConfig={this.SetAppConfig}/>

                        */}
                    </Panel>
            </div>
        );
    }

});

var WelcomeTreeTimePage = React.createClass({

    componentDidMount: function(){
        var parentNode = this.getDOMNode().parentNode;
        // set user id from the HTML template
        this.setAppState({'UID':parentNode.attributes.userid.value});
    },

    setAppState : function(partialState){
        this.setState(partialState);
    },

    render:function(){
        return (
            <div>
                <Header/>
                <PanelText/>
                <div className='bigspacer'/>
                <PanelFiles/>

                <PanelExamples/>

                <PanelConfig/>
                <Button bsStyle="primary">Run treetime</Button>
            </div>
        );
    }
});

ReactDOM.render((
    <WelcomeTreeTimePage/>),
    document.getElementById('react'));

export default WelcomeTreeTimePage;
