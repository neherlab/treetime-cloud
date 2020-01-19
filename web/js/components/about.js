import React from  'react'
import ReactDOM from 'react-dom'
import { Panel, Button, Grid, Row, Col, Glyphicon } from "react-bootstrap";

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
                                Department
                            </Col>

                            <Col xs={10} md={8} style={{"text-align":"left"}}>
                                Biozentrum
                            </Col>

                             <Col xs={6} md={4} style={{"text-align":"right"}}>
                                Institution
                            </Col>

                            <Col xs={10} md={8} style={{"text-align":"left"}}>
                                University of Basel
                            </Col>

                            <Col xs={6} md={4} style={{"text-align":"right"}}>
                                Address
                            </Col>

                            <Col xs={10} md={8} style={{"text-align":"left"}}>
                                Klingelbergstrasse 50/70, 4056 Basel, Switzerland
                            </Col>
                            <Col xs={6} md={4} style={{"text-align":"right"}}>
                                Email
                            </Col>

                            <Col xs={10} md={8} style={{"text-align":"left"}}>
                                richard.neher _AT_ unibas.ch
                            </Col>
                            <Col xs={6} md={4} style={{"text-align":"right"}}>
                                Website
                            </Col>

                            <Col xs={10} md={8} style={{"text-align":"left"}}>
                                <a href="https://neherlab.org/richard-neher.html" target="_blank" title="Personal website of Richard Neher"><Glyphicon glyph="link"/></a>
                            </Col>

                        </Row>
                    </Col>
                </Row>
                <div className="bigspacer"/>
                <Row style={{"margin-top":"40px", "margin-bottom":"40px"}}>
                    <Col xs={4} md={2}>
                    <img src="/static/svg/pavsag.jpg" width="150" />
                    </Col>

                    <Col xs={12} md={8}>
                        <Row>
                            <h3 className="text-primary" style={{"text-align":"center"}}>Pavel Sagulenko</h3>
                        </Row>

                        <Row>

                             <Col xs={6} md={4} style={{"text-align":"right"}}>
                                Department
                            </Col>

                            <Col xs={10} md={8} style={{"text-align":"left"}}>
                                Evolutionary Dynamics and Biophysics group
                            </Col>

                             <Col xs={6} md={4} style={{"text-align":"right"}}>
                                Institution
                            </Col>

                            <Col xs={10} md={8} style={{"text-align":"left"}}>
                                Max Planck Institute for Developmental Biology
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
                                <a href="http://neherlab.org/pavel-sagulenko.html" target="_blanck" title="Personal website of Pavel Sagulenko"><Glyphicon glyph={"link"}/></a>
                            </Col>

                        </Row>
                    </Col>
                </Row>
                <Row>

                <Col  xs={3} md={4}>
                <a href="http://erc.europa.eu/intra-patient-evolution-hiv" target="_blank" title="European Research Council">
                            <img src="/static/svg/logo-erc.jpg" width="100" /></a>
                </Col>

                <Col  xs={3} md={3}>
                 <a href="http://www.mpg.de" target="_blank" title="Max-Planck-Gesellschaft">
                            <img src="/static/svg/mpg.png" width="120" /></a>
                </Col>

                <Col  xs={3} md={4}>
                 <a href="http://www.unibas.ch" target="_blank" title="University of Basel">
                            <img src="/static/svg/unibas_logo.png" width="180" /></a>
                </Col>

                </Row>

                </Grid>
            </Panel>

            <Panel collapsible defaultCollapsed header="Source code and issues">
                <div style={{"text-align":"left"}}>
                    The treetime code is available under an MIT license on <a href="https://github.com/neherlab/treetime">Github</a>. Please submit any issues there.
                </div>
                <div style={{"text-align":"left"}}>
                    The web server is built using <a href="http://flask.pocoo.org/" target="_blank">flask</a>. Our code can be found at <a href="https://github.com/neherlab/treetime_web">github.com/neherlab/treetime_web</a>.
                </div>
            </Panel>

            <Panel collapsible defaultCollapsed header="Haftung & Copyright">
            <div  style={{"text-align":"justify"}}>
<h4>Web Analytics</h4> This website uses google analytics to collect statistics on
page usage. To this end, cookies are stored on your computer and javascript is
executed by your browser. If you object to such tracking tools, you can block
them using appropriate browser plug-ins.


<h4>Haftungsausschluss</h4>

Das Biozentrum der Universität Basel bemüht sich um richtige Informationen auf
seiner Homepage. Die Universität Basel übernimmt keinerlei Gewähr oder
Zusicherungen für die Vollständigkeit, Zuverlässigkeit, Aktualität,
Genauigkeit und inhaltliche Richtigkeit der bereitgestellten Informationen.

Der Zugang, die Benutzung der Homepage und die Inanspruchnahme der
Dienstleistungen geschehen auf eigenes Risiko des Benutzers / der Benutzerin.
Weder die Universität Basel noch eine von ihr beigezogene Hilfsperson, die bei
der Herstellung, Informationseingabe und -weitergabe dieser Homepage
involviert sind, sind haftbar für irgendwelche Schäden materieller oder
immaterieller Art und Investitionen, die aus dem Zugang, der Benutzung bzw.
Nichtnutzung der veröffentlichten Inhalte oder allfälligen Störungen im
Gebrauch der Homepage entstanden sind.

Das Biozentrum der Universität Basel behält sich ausdrücklich das Recht vor,
jederzeit aus welchem Grund und in welcher Art und Weise auch immer den Inhalt
dieser Homepage ohne vorgehende Ankündigung zu ändern, zu löschen oder
zeitweise nicht zu veröffentlichen. Das Biozentrum der Universität Basel lehnt
jegliche Verantwortung für irgendwelche Folgen aus solchen Änderungen,
Löschungen oder Nichtveröffentlichungen ab.

<h4>Urheberrechte</h4>

Sämtliche Online-Inhalte (Dokumente, Webseiten und deren Teile) auf der
Homepage des Biozentrums der Universität Basel sind urheberrechtlich geschützt
und dürfen nur zum privaten, wissenschaftlichen und nicht kommerziellen
Gebrauch kopiert und ausgedruckt werden.

Jegliche Vervielfältigung, Wiedergabe, Weiterübertragung oder sonstiger
Gebrauch dieser Online-Inhalte für kommerzielle Zwecke ist untersagt. Dies
gilt insbesondere auch für das Logo der Universität Basel. Allfällige
Bewilligungsgesuche sind schriftlich an die Abteilung Kommunikation des
Biozentrums der Universität Basel zu richten.

<h4>Liability</h4>

The Biozentrum of the University of Basel makes every effort to ensure that
the information on its homepage is correct. The University of Basel gives no
guarantee for the completeness, reliability, topicality, accuracy or
correctness of the information provided.

Users access and make use of the website, and any related services, at their
own risk. Neither the University of Basel nor any auxiliary person recruited
by the university to help with the production, entering and release of
information used on this website can be held liable for any material or non-
material damage, or for expenditure incurred in connection with access to, or
use or non-use of published content, or possible interruptions in the
operation of the website.

The Biozentrum at the University of Basel specifically reserves the right to
amend, delete or not publish website content for whatever reason, and in
whatever manner, without prior notification. The Biozentrum at the University
of Basel cannot accept any responsibility for repercussions of any kind
arising out of such amendments, deletions or non-publication.

The University of Basel does not check third-party websites, i.e. websites
which are not hosted on its servers, or lie beyond its scope of influence, but
which may be connected to the homepage of the Biozentrum of the University of
Basel via hyperlinks; the University of Basel declines all responsibility for
the content of such websites, or for products or services offered via such
websites.

<h4>Copyright</h4>

All online content (documents, websites and parts thereof) appearing on the
website of the Biozentrum of the University of Basel is protected by copyright
and may only be copied and/or printed out for private, scientific, non-commercial use.

All forms of reproduction, playback, onward transfer or other use of this
online content for commercial purposes is prohibited. This also applies to use
of the University of Basel logo in particular. Please contact the Biozentrum
Communications Department at the University of Basel to apply in writing for
permission to use online content.

Please note: Only the German version of the Impressum/Disclaimer is legally
binding (see above). The English translation thereof is provided for
convenience purposes only and is not legally binding.
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
