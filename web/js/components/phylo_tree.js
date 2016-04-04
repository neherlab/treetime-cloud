var EventEmitter = require('events').EventEmitter;
var d3 = require('d3');

var Globals = require('./globals.js')
var colors = Globals.colors;
var ANIMATION_DURATION = Globals.ANIMATION_DURATION;

var PhyloTree = {};

PhyloTree._lastState = null;

PhyloTree.create = function(el, props, state){

    console.log("CREATING tree");

    var svg = d3.select(el).append('svg')
      .attr('class', 'd3')
      .attr('width', props.width)
      .attr('height', props.height);
    
    svg.append('g')
        .attr('class', 'd3-year-tick')

    svg.append('g')
        .attr('class', 'd3-year-grid')
    

    svg.append('g')
        .attr('class', 'd3-links');

    svg.append('g')
        .attr('class', 'd3-tips');

    var dispatcher = new EventEmitter();
    this._create_data(props);
    
    this._update_node_colors(state);
//    this.update(el, state, dispatcher);
    return dispatcher;

};

PhyloTree._create_data = function(props){
    
    console.log(props.root)
    this.root = props.root;
    this.layout = d3.layout.tree();
    this.all_nodes = this.layout.nodes(this.root);
    this.all_links = this.layout.links(this.all_nodes);
    this.all_links.push({"source":this.root, "target":this.root});
    this.all_tips = this.gatherTips(this.root, []);
    this.all_tips.map(function(d){d.selected=false;});
    // axis
    this.year_ticks = [];
    this.month_ticks = [];

};

PhyloTree.gatherTips = function (node, tips) {
    
    if (typeof node.children != "undefined") {
        for (var i=0, c=node.children.length; i<c; i++) {
            this.gatherTips(node.children[i], tips);
        }
    }
    else {
        tips.push(node);
    }
    return tips;

};

PhyloTree.update = function(el, state, dispatcher) {
  
  console.log("UPDATING tree");

  this._update_node_colors(state);
  
  this.all_tips.map(function(d){
      if (state.selected_tip && d.strain == state.selected_tip){
          d.selected = true;
      }else{
          d.selected = false;
      }
  });
  
  this._update_scales(el, state, dispatcher);

};

PhyloTree._update_scales = function(el, state, dispatcher){
    var scales = this._scales(el, state);
    this._update_axis (el, state, scales, dispatcher);
    this._set_node_coordinates(scales, state);
    this._drawLinks(el, scales, dispatcher);
    this._drawTips(el, scales, dispatcher); 
};

PhyloTree._scales = function(el, state) {

  var width = el.offsetWidth;
  var height = el.offsetHeight;
  console.log("W = " + width);
  
  var padding = 10;

  var xUnit = state.xUnit  
  var xValues = this.all_nodes.map(function(d) {return +d[xUnit];});
  var yValues = this.all_nodes.map(function(d) {return +d.yvalue;});
  console.log(d3.min(xValues))
  var x = d3.scale.linear()
      .domain([d3.min(xValues), d3.max(xValues)])
      .range([padding, width-padding]);
  
  var y = d3.scale.linear()
      .domain([d3.min(yValues), d3.max(yValues)])
      .range([padding,height-padding])

  return {x: x, y: y};

};

PhyloTree._update_node_colors = function(state){
    
    this.all_nodes.forEach(function (d) {
        d.color = state.cscale.get_color(d[state.cscale.cUnit]);
    });

};

PhyloTree.get_color_scale = function(){

};

PhyloTree._tipFillColor = function(d) {
    
    var c =  d3.rgb(d.color).brighter([0.65]);
    return c;
}

PhyloTree._tipStrokeColor = function(d)  {return d.color;}

PhyloTree._branchStrokeColor = function (d) {return "#BBBBBB";}

PhyloTree._tipVisibility = function (d) { return d.current?"visible":"hidden";}

PhyloTree._branchStrokeWidth = function(d) {return 3;}

PhyloTree._set_node_coordinates = function(scales, state){

    var xUnit = state.xUnit
    this.all_nodes.forEach(function (d) {
            d.old_x = d.x; // need this to make smooth transitions 
            d.old_y = d.y;
            
            d.x = scales.x(d[xUnit]);
            d.y = scales.y(d.yvalue);
    });

};

PhyloTree._tipRadius = function(d) {return d.selected ? 6.0 : 4.0;}

PhyloTree._drawTips = function(el, scales, dispatcher) {
    
    var g = d3.select(el).selectAll('.d3-tips');

    var tip = g.selectAll('.d3-tip')
        .data(this.all_tips);

    tip.enter()
      .append("circle")
      .attr("class", "d3-tip")
      .attr("id", function(d) { return (d.strain).replace(/\//g, ""); })
      .attr("cx", function(d) { 
            if (d.old_x){
                return d.old_x;    
            }
            return d.x; 
        })

    tip
      .attr("r", this._tipRadius)
      .attr("cy", function(d) {return d.y;})
      .style("fill", this._tipFillColor)
      .style("stroke", this._tipStrokeColor)
      .on('mouseover', function(d) {
          console.log("Mouseover" + d.y);
          dispatcher.emit('point:tip_mouseover', d);
      })
      .on('mouseout', function(d) {
          dispatcher.emit('point:tip_mouseout', d);
      })
      .transition()
      .duration(ANIMATION_DURATION)
      .attr('cx', function(d) { return d.x; })
      .attr('cy', function(d){return d.y; });

      //.style("stroke", this._tipStrokeColor);

};

PhyloTree._drawLinks = function(el, scales, dispatcher){
        
        var g = d3.select(el).selectAll('.d3-links');

        var link = g.selectAll('.d3-link')
            .data(this.all_links);

        link.enter()
            .append("polyline")
            .attr("class", "d3-link")
            .attr("points", this._branchPoints_old)
            .style("stroke-width", 3)
            .style("stroke", "gray")
            .style("stroke-linejoin", "round")
            .style("cursor", "pointer");
            
        link.transition()
                .attr("points", this._branchPoints)
                .duration(ANIMATION_DURATION)
            
};

PhyloTree._branchPoints = function(d) {
    var tmp =   d.source.x.toString() + "," + d.target.y.toString() + " "
              + d.target.x.toString() + "," + d.target.y.toString();
    if (typeof d.target.children != "undefined"){
        var child_ys = d.target.children.map(function (x){return x.y;});
        tmp+= " "+ d.target.x.toString()+","+d3.min(child_ys).toString() + " "
                 + d.target.x.toString()+","+d3.max(child_ys).toString();
    }
    return tmp;
},

PhyloTree._branchPoints_old = function(d) {
    var tmp =   d.source.old_x.toString() + "," + d.target.old_y.toString() + " "
              + d.target.old_x.toString() + "," + d.target.old_y.toString();
    if (typeof d.target.children != "undefined"){
        var child_ys = d.target.children.map(function (x){return x.old_y;});
        tmp+= " "+ d.target.old_x.toString()+","+d3.min(child_ys).toString() + " "
                 + d.target.old_x.toString()+","+d3.max(child_ys).toString();
    }
    return tmp;
},


PhyloTree._update_axis = function(el, state, scales, dispatcher){

    if (state.xUnit == "xvalue"){
        this._hide_axis(el, state, scales, dispatcher);
    }else{
        this._draw_axis(el, state, scales, dispatcher);
    }

};

PhyloTree._draw_axis = function(el, state, scales, dispatcher){
    
    var g = d3.select(el).selectAll('.d3-year-grid');
    
    var bnd = scales.x.domain()
    var yticks_xs = d3.range(bnd[0], bnd[1], 1);
    var y_bnds = scales.y.range()
    var glines = yticks_xs.map(function(d){return ({x1:d, y1:y_bnds[0], x2:d, y2:y_bnds[1]});})

    console.log(glines)
    var gridLine = g.selectAll('.d3-year-grid-line')
        .data(glines);

    gridLine.enter()
        .append('line')
        .attr("class","d3-year-grid-line")
        .attr("x1", function(d){return scales.x(d.x1);})
        .attr("x2", function(d){return scales.x(d.x2);})
        .attr("y1", function(d){return d.y1;})
        .attr("y2", function(d){return d.y2;})
        .style("stroke-width", 1)
        .style("stroke", "#DDDDDD");
};

PhyloTree._hide_axis = function (el, state, scales, dispatcher){
    var g = d3.select(el).selectAll('.d3-year-grid-line');
    g.remove();

};

PhyloTree.destroy = function(el) {

};

export default PhyloTree;
