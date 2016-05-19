import React from  'react'
import ReactDOM from 'react-dom'
import { Panel, Button, Grid, Row, Col } from "react-bootstrap";

import Header from './header.js'
import Footer from './footer.js'

var Terms = React.createClass({
    render: function(){
        return (
            <div>
            <Header/>
            <div className="page_container">
            .
            <Panel collapsible defaultExpanded header="People and funding">
                <Grid> 
                
                
                <Row> 
                    <Col xs={4} md={2}> 
                    <img src="/static/svg/rneher.jpg" width="150" />
                    </Col> 
                    
                    <Col xs={12} md={8}>
                        <Row> 
                            <h3 className="text-primary">Richard Neher</h3>
                        </Row> 

                        <Row>
                            <Col xs={6} md={4} style={{"text-align":"right"}}> 
                                Institute
                            </Col> 

                            <Col xs={10} md={8} style={{"text-align":"left"}}> 
                                Max Planck Institute
                            </Col> 

                             <Col xs={6} md={4} style={{"text-align":"right"}}> 
                                Department
                            </Col> 

                            <Col xs={10} md={8} style={{"text-align":"left"}}> 
                                Developmental Biology
                            </Col> 

                            <Col xs={6} md={4} style={{"text-align":"right"}}> 
                                Address
                            </Col> 

                            <Col xs={10} md={8} style={{"text-align":"left"}}> 
                                Spemannstrasse 35, 72076 Tuebingen, Germany
                            </Col> 
                            <Col xs={6} md={4} style={{"text-align":"right"}}> 
                                Email
                            </Col> 

                            <Col xs={10} md={8} style={{"text-align":"left"}}> 
                                richard.neher _AT_ tuebingen.mpg.de
                            </Col>
                            <Col xs={6} md={4} style={{"text-align":"right"}}> 
                                Website
                            </Col> 

                            <Col xs={10} md={8} style={{"text-align":"left"}}> 
                                <a href="http://www.eb.tuebingen.mpg.de/nc/research/departments/details/details/rneher.html" target="_blanck" title="Personal website of Richard Neher"><span class="glyphicon glyphicon-link"></span></a>
                            </Col>

                        </Row> 
                    </Col> 
                </Row>

                <Row> 
                    <Col xs={4} md={2}> 
                    <img src="/static/svg/pausag.jpg" width="150" />
                    </Col> 
                    
                    <Col xs={12} md={8}>
                        <Row> 
                            <h3 className="text-primary" style={{"text-align":"center"}}>Pavel Sagulenko</h3>
                        </Row> 

                        <Row>
                            <Col xs={6} md={4} style={{"text-align":"right"}}> 
                                Institute
                            </Col> 

                            <Col xs={10} md={8}  style={{"text-align":"left"}}> 
                                Max Planck Institute
                            </Col> 

                             <Col xs={6} md={4} style={{"text-align":"right"}}> 
                                Department
                            </Col> 

                            <Col xs={10} md={8} style={{"text-align":"left"}}> 
                                Developmental Biology
                            </Col> 

                            <Col xs={6} md={4} style={{"text-align":"right"}}> 
                                Address
                            </Col> 

                            <Col xs={10} md={8} style={{"text-align":"left"}}> 
                                Spemannstrasse 35, 72076 Tuebingen, Germany
                            </Col> 
                            <Col xs={6} md={4} style={{"text-align":"right"}}> 
                                Email
                            </Col> 

                            <Col xs={10} md={8} style={{"text-align":"left"}}> 
                                pavel.sagulenko _AT_ tuebingen.mpg.de
                            </Col>
                            <Col xs={6} md={4} style={{"text-align":"right"}}> 
                                Website
                            </Col> 

                            <Col xs={10} md={8} style={{"text-align":"left"}}> 
                                <a href="http://www.eb.tuebingen.mpg.de/nc/research/departments/details/details/rneher.html" target="_blanck" title="Personal website of Richard Neher"><span class="glyphicon glyphicon-link"></span></a>
                            </Col>

                        </Row> 
                    </Col> 
                </Row>
                <Row>
                
                <Col  xs={6} md={4}> 
                <a href="http://erc.europa.eu/intra-patient-evolution-hiv" target="_blank" title="European Research Council">
                            <img src="/static/svg/logo-erc.jpg" width="100" /></a>
                </Col> 
                
                <Col  xs={6} md={4}> 
                 <a href="http://www.mpg.de" target="_blank" title="Max-Planck-Gesellschaft">
                            <img src="/static/svg/mpg.png" width="120" /></a>
                </Col> 
                
                </Row>

                </Grid> 
            </Panel>      
            
            <Panel collapsible defaultExpanded header="Source code and issues">
                <div style={{"text-align":"left"}}> 
                    Please visit the <a href="https://github.com/neherlab/treetime">Github</a> project, which is open source.
                    </div> 
            </Panel>            
            
            <Panel collapsible defaultCollapsed header="Data protection advice">
            <div  style={{"text-align":"justify"}}>
            <h4 class="text-muted">Data Collection and Processing</h4>
          
          We wish to make our websites as appealing and comfortable for you as possible. To this end, analysis of statistical information relating to your utilisation of our websites and the collection of technical details on the web browsers and computers of our online visitors is very helpful for us.

          We employ the web analytics tool Piwik for these purposes, which uses Cookies and JavaScript to collect information on your computer and automatically forwards it to us (pseudonymised user data).

          Storing and analysis of data takes place on a Piwik server which is operated by the Max Planck Institute for Developmental Biology. Of course, you may object to the collection of data. Objecting is very simple, it only takes one click. Further information on which data is collected and on your options to object is available in our Data Collection Declaration at the end of this page.

          Also, the integration of external services such as Google Maps for route maps, Youtube or Amazon Cloud is always undertaken in a considerate manner and with the aim of making your visit on our websites as pleasant as possible.

          We must advise you, however, that your IP address and perhaps other data related to your person will be transmitted to the service provider concerned (Google, Amazon etc.) and may be stored or analyzed there.

          The data protection officer of the Max Planck Society is available for questions concerning the topic of data protection at dsb[at]gv.mpg.de or +49 89 2108-1554.


          <h4 class="text-muted">Data Transmission</h4>
          Your personal data will only be transmitted to government organizations and authorities in legally required cases and/or for prosecution in the event of attacks on our network infrastructure. Your personal data are not provided to third parties for any other purpose.

          <h4 class="text-muted">Liability for Contents of Online Information</h4>
          As the provider of contents in accordance with Section 7 Paragraph 1 of the Tele-Media Law, the Max Planck Society shall be responsible for any contents which it makes available for use in accordance with general legal provisions. The Max Planck Society makes every effort to provide timely and accurate information on this Web site. Nevertheless, errors and inaccuracies cannot be completely ruled out. Therefore, the Max Planck Society does not assume any liability for the relevance, accuracy, completeness or quality of the information provided. The Max Planck Society shall not be liable for damage of a tangible or intangible nature caused directly or indirectly through the use or failure to use the information offered and/or through the use of faulty or incomplete information unless it is verifiably culpable of intent or gross negligence. The same shall apply to any downloadable software available free of charge. The Max Planck Society reserves the right to modify, supplement, or delete any or all of the information offered on its Internet site, or to temporarily or permanently cease publication thereof without prior and separate notification.  


          <h4 class="text-muted">Links to Internet Sites of Third Parties</h4>
          This Internet site includes links to external pages. These external links are designated as follows: [blue bullet, underlined, green]. References to the subordinate pages of this Internet site are designated as follows: [straight arrow ].  
          The respective provider shall be responsible for the contents of any linked external pages. In establishing the initial link, the Max Planck Society has reviewed the respective external content in order to determine whether such link entailed possible civil or criminal responsibility. However, a constant review of linked external pages is unreasonable without concrete reason to believe that a violation of the law may be involved. If the Max Planck Society determines such or it is pointed out by others that an external offer to which it is connected via a link entails civil or criminal responsibility, then the Max Planck Society will immediately eliminate any link to this offer. The Max Planck Society expressly dissociates itself from such contents.  

          <h4 class="text-muted">Copyright</h4>
         
          The layout, graphics employed and any other contents on the homepage of the Max Planck Society Internet site are protected by copyright law. © Max-Planck-Gesellschaft zur Förderung der Wissenschaften e.V., Munich. All rights reserved.  
         </div>
            </Panel>
            </div>
            <Footer/>
            </div>
            );
    }
});

ReactDOM.render((
    <Terms />),
    document.getElementById('react'));


export default Terms;